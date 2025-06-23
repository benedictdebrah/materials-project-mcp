# Materials Project MCP

A Model Context Protocol (MCP) server for querying the Materials Project database using the mp_api client.

## Requirements

- **Materials Project API Key** - [Get one here](https://materialsproject.org/) (free account required)
- **Docker Desktop** (must be running) OR **Python 3.12+** with [uv](https://github.com/astral-sh/uv)



### Getting Your Materials Project API Key

1. Visit [Materials Project](https://materialsproject.org/)
2. Create a free account or log in
3. Go to your dashboard
4. Navigate to API settings
5. Generate or copy your API key
6. Keep this key secure - you'll need it for setup

## Installation Options

### Step 1: Docker (Recommended)

#### Using Docker Run
1. **Install Docker Desktop:**
   - Download from [docker.com](https://www.docker.com/products/docker-desktop/)
   - Install and **make sure Docker Desktop is running**

2. **Pull the Docker image:**
   ```bash
   docker pull benedict2002/materials-project-mcp
   ```

3. **Test the installation:**
   ```bash
   docker run --rm -i -e MP_API_KEY="your-api-key" benedict2002/materials-project-mcp
   ```

#### Using Docker Compose (Easiest)
1. **Install Docker Desktop and make sure it's running**

2. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd materials-project-mcp
   ```

3. **Create a `.env` file:**
   ```bash
   echo "MP_API_KEY=your-materials-project-api-key" > .env
   ```

4. **Test the setup:**
   ```bash
   docker-compose up
   ```

5. **For background running:**
   ```bash
   docker-compose up -d
   ```

6. **Stop the service:**
   ```bash
   docker-compose down
   ```

### Step 1 : Local Python Installation

1. **Install uv (if not already installed):**
   ```bash
   curl -Ls https://astral.sh/uv/install.sh | sh
   ```

2. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd materials-project-mcp
   ```

3. **Create and activate virtual environment:**
   ```bash
   uv venv
   source .venv/bin/activate  # Linux/macOS
   # or
   .venv\Scripts\activate     # Windows
   ```

4. **Install dependencies:**
   ```bash
   uv pip install -r requirements.txt
   ```

5. **Set your API key:**
   ```bash
   export MP_API_KEY="your-api-key"  # Linux/macOS
   # or
   set MP_API_KEY=your-api-key       # Windows
   ```

6. **Test the installation:**
   ```bash
   python server.py
   ```

## Step 2 : Setup with Claude Desktop

1. **Locate your Claude configuration file:**
   - **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

2. **Choose your configuration method:**

   **Using Docker Run**
   ```json
   {
     "mcpServers": {
       "Materials Project MCP": {
         "command": "docker",
         "args": [
           "run", "--rm", "-i",
           "-e", "MP_API_KEY=your-materials-project-api-key",
           "benedict2002/materials-project-mcp"
         ]
       }
     }
   }
   ```

3. **Replace `your-materials-project-api-key` with your actual API key**

4. **Ensure Docker Desktop is running**

5. **Restart Claude Desktop**

6. **Verify installation:**
   - Open a new chat in Claude
   - Ask something like "Search for silicon materials in the Materials Project database" or test any of the availabe tools.
   - You should see Materials Project data in the response

## Setup with VS Code Copilot

1. **Open VS Code Settings:**
   - Press `Ctrl+Shift+P` (Windows/Linux) or `Cmd+Shift+P` (macOS)
   - Type "Preferences: Open User Settings (JSON)"
   - Select it to open settings.json

2. **Add MCP configuration:**
   ```json
   {
     "mcp": {
       "inputs": [],
       "servers": {
         "Materials Project MCP": {
           "command": "docker",
           "args": [
             "run", "--rm", "-i",
             "-e", "MP_API_KEY=your-api-key",
             "benedict2002/materials-project-mcp"
           ]
         }
       }
     },
     "chat.mcp.discovery.enabled": true,
     "workbench.secondarySideBar.showLabels": false
   }
   ```

3. **Alternative: Local Python setup for VS Code:**
   ```json
   {
     "mcp": {
       "inputs": [],
       "servers": {
         "Materials Project MCP": {
           "command": "/usr/local/bin/uv",
           "args": [
             "run",
             "--with",
             "mcp[cli],aiohttp,pydantic,mp_api,pymatgen,emmet-core",
             "/path/to/your/server.py"
           ],
           "env": {
             "MP_API_KEY": "your-api-key"
           }
         }
       }
     },
     "chat.mcp.discovery.enabled": true
   }
   ```

4. **Replace placeholders:**
   - `your-api-key` with your Materials Project API key
   - `/path/to/your/server.py` with the actual path to server.py

5. **Ensure Docker Desktop is running** (for Docker configurations)

6. **Restart VS Code**

7. **Test in VS Code:**
   - Open VS Code chat/copilot
   - Ask about materials from the Materials Project
   - The Docker container will start automatically when VS Code makes requests


## Testing & Development (developers)

### Testing Your Installation

1. **Test MCP server locally:**
   ```bash
   mcp dev server.py
   ```
   Look for the line "🔗 Open inspector with token pre-filled:" and use that URL

### Development Workflow

1. **Create a feature branch:**
   ```bash
   git checkout -b feature-name
   ```

2. **Make your changes and test:**
   ```bash
   # Local testing with MCP Inspector
   mcp dev server.py
   # Use the inspector URL to test your changes interactively
   
   # Docker testing
   docker build -t materials-project-mcp-local .
   docker run --rm -i -e MP_API_KEY="your-api-key" materials-project-mcp-local
   
   # Docker Compose testing
   docker-compose up --build
   ```

3. **Commit and push:**
   ```bash
   git add .
   git commit -m "Add feature description"
   git push origin feature-name
   ```

4. **Open a pull request**


## Available Tools & Features

- **search_materials** - Search by elements, band gap range, stability
- **get_structure_by_id** - Get crystal structures and lattice parameters
- **get_electronic_bandstructure** - Plot electronic band structures
- **get_electronic_dos_by_id** - Get electronic density of states
- **get_phonon_bandstructure** - Plot phonon band structures
- **get_phonon_dos_by_id** - Get phonon density of states
- **get_ion_reference_data_for_chemsys** - Download aqueous ion reference data for Pourbaix diagrams
- **get_cohesive_energy** - Calculate cohesive energies
- **get_atom_reference_data** - Retrieve reference energies of isolated neutral atoms
- **get_magnetic_data_by_id** - Magnetic properties and ordering
- **get_charge_density_by_id** - Charge density data
- **get_dielectric_data_by_id** - Dielectric constants and properties
- **get_diffraction_patterns** - X-ray and neutron diffraction
- **get_xRay_absorption_spectra** - XAFS, XANES, EXAFS spectra
- **get_elastic_constants** - Mechanical properties
- **get_suggested_substrates** - Find substrates for thin films
- **get_thermo_stability** - Thermodynamic stability analysis
- **get_surface_properties** - Surface energies, work functions, and Wulff shapes
- **get_grain_boundaries** - Computed grain boundaries for a material
- **get_insertion_electrodes** - Insertion electrode and battery data
- **get_oxidation_states** - Element oxidation states, formula, and structure info

## Troubleshooting

### Common Issues

1. **"Invalid API key" error:**
   - Verify your API key is correct
   - Check that you've set the environment variable properly
   - Ensure your Materials Project account is active

2. **"Docker not found" or "Cannot connect to Docker daemon":**
   - **Make sure Docker Desktop is installed and running**
   - You should see the Docker Desktop icon in your system tray/menu bar
   - Try `docker --version` to verify Docker is accessible
   - On Windows/Mac: Open Docker Desktop application
   - On Linux: Start Docker service with `sudo systemctl start docker`

3. **Container startup issues:**
   - Docker containers start automatically when Claude/VS Code makes requests
   - No need to manually start containers - they're ephemeral (start → run → stop)
   - Each query creates a fresh container instance

4. **Docker Compose issues:**
   - Make sure Docker Compose is installed: `docker-compose --version`
   - Check your `.env` file exists and has the correct API key
   - Verify the docker-compose.yml file is in the correct location
   - Ensure Docker Desktop is running

5. **MCP server not recognized in Claude:**
   - Check your configuration file path
   - Verify JSON syntax is correct
   - Restart Claude Desktop after configuration changes
   - Ensure Docker Desktop is running


### Getting Help

- **MCP Inspector**: Use `mcp dev server.py` for interactive testing and debugging.
- **GitHub Issues**: [Create an Issue](https://github.com/yourusername/materials-project-mcp/issues) for bug reports, feature requests, or questions.
- **Materials Project API Docs**: [docs.materialsproject.org](https://docs.materialsproject.org/)
- **MCP Documentation**: [modelcontextprotocol.io](https://modelcontextprotocol.io/)
- **Docker Help**: [docs.docker.com](https://docs.docker.com/)

---
### Authors

- Benedict Debrah
- Peniel Fiawornu

### Reference

Yin, Xiangyu. 2025. "Building an MCP Server for the Materials Project." March 23, 2025. https://xiangyu-yin.com/content/post_mp_mcp.html.