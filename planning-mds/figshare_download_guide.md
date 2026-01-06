# Lessons Learned: Downloading Datasets from Figshare and Figshare+

This document outlines the most effective strategies for downloading scientific datasets from Figshare, especially when working in an agentic environment via SSH or remote connections.

## 1. Avoid Simple `curl` or `wget`
Basic terminal commands often trigger **403 Forbidden** errors on Figshare. This is because:
*   **Security Measures**: Figshare expects a real browser `User-Agent`.
*   **Redirect Handling**: The initial download URL usually redirects to an AWS S3 bucket with temporary credentials. Basic tools sometimes fail to carry over the necessary headers during these redirects.
*   **Domain Sensitivity**: Figshare+ (`plus.figshare.com`) may have stricter requirements than standard Figshare.

## 2. Most Reliable Method: Python `requests`
A Python script is the most robust way to download data. It provides fine-grained control over headers, automatic redirect handling, and streaming for large files.

### Recommended Script Pattern:
```python
import os
import requests
from concurrent.futures import ThreadPoolExecutor

# Use a realistic User-Agent to avoid initial 403s
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) ..."
}

def download_file(file_info):
    url = f"https://plus.figshare.com/ndownloader/files/{file_info['id']}"
    target_path = os.path.join("./data", file_info['name'])
    
    with requests.get(url, headers=HEADERS, stream=True) as r:
        r.raise_for_status() # Automatically catch 403/404 errors
        with open(target_path, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print(f"Downloaded: {file_info['name']}")

# Use parallel threads for faster batch downloads
with ThreadPoolExecutor(max_workers=5) as executor:
    executor.map(download_file, file_list)
```

## 3. Using the Browser Subagent for Diagnostics
If direct links or Python scripts fail, use the **Browser Subagent** to:
1.  Navigate to the dataset page.
2.  Capture direct S3 URLs from the network logs or by inspecting redirects via JavaScript (`XHR` or `fetch`).
3.  Verify the existence of files and their sizes.

*Note: While the browser can trigger downloads, the resulting files may be saved in an isolated "Downloads" directory that is harder to locate than a direct save path.*

## 4. Batch Identification
Always check for a metadata JSON (e.g., `replogle_20029387.json`). These often contain the file IDs and names needed to construct download loops programmatically.

## 5. Summary of Failures to Watch For
*   **0-byte files**: Indicates a silent 403 or failed redirect that was handled as a "success" by the command line tool.
*   **HTML in expected binary file**: If a downloaded `.h5ad` or `.zip` file is small and contains HTML tags, it's actually an error page (e.g., "403 Forbidden"). Check with `head -c 100 <file>`.
