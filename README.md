# Materials Project MCP

A Model Context Protocol (MCP) server that exposes query tools for the Materials Project database using the mp_api client.


## Getting Started

1. **Clone this repository:**


2. **Make sure you have [uv](https://github.com/astral-sh/uv) installed:**

   ```bash
   curl -Ls [https://astral.sh/uv/install.sh](https://astral.sh/uv/install.sh) | sh
   ```

3. **Create a virtual environment using `uv`:**

   ```bash
   uv venv
   source .venv/bin/activate  # or .venv\Scripts\activate on Windows
   ```

4. **Install dependencies from `requirements.txt`:**

   ```bash
   uv pip install -r requirements.txt
   ```

---

## 🛠️ Development Workflow

### Branching

* When working on a new feature or task, **create a new branch** named after the tool or feature you're working on.

  ```bash
  git checkout -b <tool-name>
  ```

### Pushing Changes

* After completing your task and committing your changes, **push your branch**:

  ```bash
  git push origin <tool-name>
  ```

### Pull Requests

* Open a **Pull Request (PR)** from your branch to `main`.
* Use clear titles and descriptions to explain the changes and context.

---

## ✅ Example Workflow

```bash
git checkout -b image-generation
# make your changes...
git add .
git commit -m "Add image generation feature"
git push origin image-generation
# Then open a pull request

```

## Usage

### Testing the Server

To test the server locally with mcp inspector, you can use the MCP development tool: check mcp docs to know more

```bash
mcp dev server.py
```

### Using with Claude

To use this server with Claude:

1. Install the server:
```bash
mcp install server.py
```

2. Edit the generated config file with the following configuration (replace the path with your actual path):
```json
{
  "mcpServers": {
    "Materials Project": {
      "command": "/usr/local/bin/uv",
      "args": [
        "run", 
        "--with", 
        "mcp[cli],aiohttp,pydantic,mp_api,pymatgen,emmet-core", 
        "/path/to/your/materials-project-mcp/server/server.py"
      ],
      "env": {
        "MP_API_KEY": "YOUR_API_KEY" 
      }
    }
  }
}
```

Note: Make sure to replace `YOUR_API_KEY` with your actual Materials Project API key.

## Features

- Search for materials by elements, band gap range, and stability
- Retrieve crystal structures by Materials Project ID
- Integration with Claude AI assistant

## Requirements

- Python 3.12+
- Materials Project API key
- Required packages: mcp, aiohttp, pydantic, mp_api



### Citations


1. Yin, Xiangyu. 2025. “Building an MCP Server for the Materials Project.” March 23, 2025. https://xiangyu-yin.com/content/post_mp_mcp.html.