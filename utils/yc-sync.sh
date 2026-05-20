#!/usr/bin/env bash
set -euo pipefail

# yc_sync.sh

# -----------------------------
# Настройки, зафиксированные в коде
# -----------------------------
ROOT_DIR="$HOME/bioinf/0_DIPLOM"

SSH_DIR="$HOME/.ssh/yc_sync"
YC_LOGIN="hilenka2026"
YC_ORG_ID="bpfndhd51nfbl3nkh08f"

LOCAL_DIR="$ROOT_DIR/local"
REMOTE_DIR="/mnt/low"

# Файл с настройками, которые создаются командой configure
CONFIG_FILE="$ROOT_DIR/.yc_sync.conf"

# -----------------------------
# Вспомогательные функции
# -----------------------------

usage() {
  cat <<EOF
Использование:

  $0 reload
      Пересоздает SSH-конфигурацию через yc compute ssh certificate export

  $0 configure -a SERVER_ADDRESS
      Сохраняет адрес удаленного сервера

  $0 get FILE_NAME
      Скачивает файл с удаленного сервера из REMOTE_DIR в LOCAL_DIR

  $0 push FILE_NAME
      Загружает файл из LOCAL_DIR на удаленный сервер в REMOTE_DIR

Примеры:

  $0 reload
  $0 configure -a 10.10.10.10
  $0 get example.txt
  $0 push example.txt
EOF
}

die() {
  echo "Ошибка: $*" >&2
  exit 1
}

load_config() {
  if [ -f "$CONFIG_FILE" ]; then
    # shellcheck disable=SC1090
    source "$CONFIG_FILE"
  fi
}

check_ready() {
  load_config

  if [ ! -d "$SSH_DIR" ]; then
    die "SSH-конфигурация не найдена. Сначала выполните: $0 reload"
  fi

  if [ -z "${SERVER_ADDRESS:-}" ]; then
    die "Адрес сервера не настроен. Сначала выполните: $0 configure -a SERVER_ADDRESS"
  fi
}

# -----------------------------
# Команды
# -----------------------------

cmd_reload() {
  DIR="$SSH_DIR"

  if [ -d "$SSH_DIR" ]; then
    rm -Rf "$SSH_DIR"
    echo "$SSH_DIR is already exists, remove"
  fi

  mkdir -p "$SSH_DIR"
  echo "Create new directory by path $SSH_DIR"

  yc compute ssh certificate export \
    --login "$YC_LOGIN" \
    --organization-id "$YC_ORG_ID" \
    --directory "$SSH_DIR"

  echo "New ssh config is created"
}

cmd_configure() {
  local server_address=""

  while getopts ":a:" opt; do
    case "$opt" in
      a)
        server_address="$OPTARG"
        ;;
      :)
        die "Опция -$OPTARG требует аргумент"
        ;;
      \?)
        die "Неизвестная опция: -$OPTARG"
        ;;
    esac
  done

  if [ -z "$server_address" ]; then
    die "Не указан адрес сервера. Используйте: $0 configure -a SERVER_ADDRESS"
  fi

  cat > "$CONFIG_FILE" <<EOF
SERVER_ADDRESS="$server_address"
EOF

  echo "Адрес сервера сохранен: $server_address"
}

cmd_get() {
  check_ready

  local file_name="${1:-}"

  if [ -z "$file_name" ]; then
    die "Не указано имя файла. Используйте: $0 get FILE_NAME"
  fi

  mkdir -p "$LOCAL_DIR"

  scp -o IdentitiesOnly=yes \
    -i "$SSH_DIR/yc-organization-id-${YC_ORG_ID}-${YC_LOGIN}" \
    "${YC_LOGIN}@${SERVER_ADDRESS}:${REMOTE_DIR}/${file_name}" \
    "${LOCAL_DIR}/${file_name}"

  echo "Файл скачан: ${LOCAL_DIR}/${file_name}"
}

cmd_push() {
  check_ready

  local file_name="${1:-}"

  if [ -z "$file_name" ]; then
    die "Не указано имя файла. Используйте: $0 push FILE_NAME"
  fi

  if [ ! -f "${LOCAL_DIR}/${file_name}" ]; then
    die "Локальный файл не найден: ${LOCAL_DIR}/${file_name}"
  fi

  scp -o IdentitiesOnly=yes \
    -i "$SSH_DIR/yc-organization-id-${YC_ORG_ID}-${YC_LOGIN}" \
    "${LOCAL_DIR}/${file_name}" \
    "${YC_LOGIN}@${SERVER_ADDRESS}:${REMOTE_DIR}/${file_name}"

  echo "Файл загружен: ${SERVER_ADDRESS}:${REMOTE_DIR}/${file_name}"
}

# -----------------------------
# Точка входа
# -----------------------------

main() {
  local command="${1:-}"

  if [ -z "$command" ]; then
    usage
    exit 1
  fi

  shift

  case "$command" in
    reload)
      cmd_reload "$@"
      ;;
    configure)
      cmd_configure "$@"
      ;;
    get)
      cmd_get "$@"
      ;;
    push)
      cmd_push "$@"
      ;;
    -h|--help|help)
      usage
      ;;
    *)
      die "Неизвестная команда: $command"
      ;;
  esac
}

main "$@"