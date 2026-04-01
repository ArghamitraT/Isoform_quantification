#!/usr/bin/env bash

set -euo pipefail

TARGET="${1:-.}"

# Resolve path
if [ -f "$TARGET" ]; then
  RESOLVED_DIR="$(cd "$(dirname "$TARGET")" && pwd)"
  FILENAME="$(basename "$TARGET")"
elif [ -d "$TARGET" ]; then
  RESOLVED_DIR="$(cd "$TARGET" && pwd)"
  FILENAME=""
else
  echo "Error: '$TARGET' is not a valid file or directory."
  exit 1
fi

# Find a free port
PORT=""
for port in 8787 8788 8789 8790; do
  if ! lsof -i :"$port" 2>/dev/null | grep -q LISTEN; then
    PORT="$port"
    break
  fi
done

if [ -z "$PORT" ]; then
  echo "Error: No free port found in 8787-8790."
  exit 1
fi

# Start server
cd "$RESOLVED_DIR"
python -m http.server "$PORT" >/dev/null 2>&1 &
PID=$!

# Build URL
if [ -n "$FILENAME" ]; then
  URL="http://localhost:$PORT/$FILENAME"
else
  URL="http://localhost:$PORT/"
fi

echo "Server PID: $PID"
echo "Server is running on port $PORT"
echo
echo "VS Code should auto-detect it in the Ports tab — click the globe icon."
echo "Or open Command Palette (Ctrl+Shift+P) -> Simple Browser: Show"
echo "Then enter:"
echo "  $URL"
echo
echo "To stop the server later:"
echo "  kill $PID"
echo
echo "Direct URL:"
echo "$URL"