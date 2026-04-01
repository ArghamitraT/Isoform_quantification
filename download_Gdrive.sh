#!/bin/bash

# === GDrive Configuration ===
GDRIVE="gdrive_columbia"
GDRIVE_BASE="technical_work/RNA_Splicing"

# === Dynamically find base Contrastive_Learning path ===
BASE_PATH=$(pwd | sed -E 's|(.*RNA_Splicing).*|\1|')

# === Dynamic local server paths ===
SERVER_DATA_BASE="${BASE_PATH}/data/final_data"
SERVER_RESULT_BASE="${BASE_PATH}/files/results"

# === Remote GDrive paths ===
GDRIVE_DATA_PATH="${GDRIVE}:${GDRIVE_BASE}/data/final_data"
GDRIVE_RESULT_PATH="${GDRIVE}:${GDRIVE_BASE}/files/results"

# === Ask for type ===
read -p "📥 Type of download? (data/result/other): " TYPE
TYPE=$(echo "$TYPE" | tr '[:upper:]' '[:lower:]')

if [[ "$TYPE" == "data" ]]; then
    echo "📁 GDrive data folders:"
    rclone lsd "$GDRIVE_DATA_PATH"
    read -p "👉 Enter folder name to download: " ITEM_NAME
    GDRIVE_PATH="${GDRIVE_DATA_PATH}/${ITEM_NAME}"
    SERVER_PATH="${SERVER_DATA_BASE}/${ITEM_NAME}"

elif [[ "$TYPE" == "result" ]]; then
    echo "📁 GDrive result folders:"
    rclone lsd "$GDRIVE_RESULT_PATH"
    read -p "👉 Enter folder name to download: " ITEM_NAME
    GDRIVE_PATH="${GDRIVE_RESULT_PATH}/${ITEM_NAME}"
    SERVER_PATH="${SERVER_RESULT_BASE}/${ITEM_NAME}"

else
    read -p "☁️  Enter full GDrive path (e.g. gdrive_columbia:technical_work/CntrstvLrn_FuncExon/data/final_data): " GDRIVE_PATH
    read -p "📂 Enter full server path to save (e.g. /path/to/save): " SERVER_PATH
fi

# === Confirm paths ===
echo ""
echo "📥 Downloading from: $GDRIVE_PATH"
echo "📁 Saving to:        $SERVER_PATH"
read -p "✅ Proceed with download? (y/n): " CONFIRM

if [[ "$CONFIRM" != "y" ]]; then
    echo "❌ Download cancelled."
    exit 0
fi

# === Check if GDRIVE_PATH is a file or folder ===
BASE_NAME=$(basename "$GDRIVE_PATH")
DIR_NAME=$(dirname "$GDRIVE_PATH")
IS_FOLDER=$(rclone lsf --dirs-only "$DIR_NAME" | grep -Fx "$BASE_NAME/")

# === Execute Download ===
if [[ -n "$IS_FOLDER" ]]; then
    echo "📁 Detected folder — using rclone copy"
    rclone copy "$GDRIVE_PATH" "$SERVER_PATH" --progress --ignore-times
else
    echo "📄 Detected file — using rclone copyto"
    rclone copyto "$GDRIVE_PATH" "$SERVER_PATH" --progress --ignore-times
fi