---
name: view-html
description: Start an HTTP server on the cluster and open an HTML file via VS Code's built-in Simple Browser. Pass a path to an HTML file or a directory.
allowed-tools: Bash
---

The user wants to view an HTML file from the cluster in their browser.

**Steps:**

1. If `$ARGUMENTS` is provided, resolve the path. If it's a file, use its parent directory as the server root and note the filename. If it's a directory, use it as the server root.

2. Find a free port (try 8787, then 8788, 8789 if taken):
```bash
for port in 8787 8788 8789 8790; do
  if ! lsof -i :$port 2>/dev/null | grep -q LISTEN; then
    echo $port; break
  fi
done
```

3. Start Python HTTP server in the background from that directory:
```bash
cd <resolved_dir> && python -m http.server <port> &
echo "Server PID: $!"
```

4. Tell the user:
   - The server is running on port `<port>`
   - VS Code should auto-detect it in the **Ports** tab (bottom panel) — click the globe icon
   - Or: open Command Palette (`Ctrl+Shift+P`) → `Simple Browser: Show` → enter `http://localhost:<port>/<filename>`
   - To stop the server later: `kill <PID>`

5. Print the direct URL: `http://localhost:<port>/<filename>`
