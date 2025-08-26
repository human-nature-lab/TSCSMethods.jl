#!/bin/bash

# Script to generate additional SVG diagrams for the visual guide
# Requires: npm install -g @mermaid-js/mermaid-cli

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSETS_DIR="$SCRIPT_DIR/src/assets/images"

# Create assets directory if it doesn't exist
mkdir -p "$ASSETS_DIR"

# Create temporary directory for mermaid files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Generating additional diagrams..."

# Advanced Workflow Diagram
cat > "$TEMP_DIR/advanced_workflow.mmd" << 'EOF'
flowchart TD
    A["✅ Basic Workflow Complete"] --> B{🚀 Advanced Options}
    
    B --> C["👥 Stratified Analysis"]
    B --> D["🎯 Refinement & Calipers"]
    B --> E["🤖 Auto-balancing"]
    B --> F["📊 Multiple Outcomes"]
    
    C --> G["stratify(model)<br/>Group-specific effects"]
    C --> H["CICStratified models<br/>Heterogeneous treatment"]
    
    D --> I["refine(model)<br/>Iterative improvement"]
    D --> J["caliper(model)<br/>Distance constraints"]
    
    E --> K["autobalance(model)<br/>P-value optimization"]
    
    F --> L["Vector outcomes<br/>Simultaneous estimation"]
    
    G --> M["📈 Subgroup Analysis<br/>Detailed Results"]
    H --> M
    I --> N["🎯 Improved Matching<br/>Better Balance"]
    J --> N
    K --> O["⚖️ Optimized Balance<br/>Statistical Precision"]
    L --> P["📊 Multi-outcome Results<br/>Joint Inference"]
    
    style B fill:#e1f5fe,stroke:#1565c0,stroke-width:2px
    style M fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style N fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style O fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style P fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
EOF

# Data Validation Diagram
cat > "$TEMP_DIR/data_validation.mmd" << 'EOF'
flowchart LR
    A["📋 Input Data<br/>DataFrame"] --> B{🔍 Data Structure Check}
    
    B --> C["⏰ Time Variable<br/>:day, :t, etc."]
    B --> D["🆔 Unit ID<br/>:fips, :id, etc."]
    B --> E["💊 Treatment<br/>:treatment, :gub"]
    B --> F["📈 Outcome<br/>:outcome, :Y"]
    B --> G["📊 Covariates<br/>:X1, :X2, :pop_dens"]
    
    C --> H["✅ Validation Checks"]
    D --> H
    E --> H
    F --> H
    G --> H
    
    H --> I["⏱️ Time Window Validation<br/>F > 0, L < 0"]
    H --> J["❓ Missing Data Handling<br/>Imputation strategies"]
    H --> K["🎯 Treatment Event Detection<br/>Binary indicators"]
    
    I --> L{All Valid?}
    J --> L
    K --> L
    
    L -->|✅ Yes| M["🚀 Proceed to Analysis<br/>Ready for matching"]
    L -->|❌ No| N["⚠️ Informative Errors<br/>Clear guidance"]
    
    style A fill:#e3f2fd,stroke:#1565c0,stroke-width:2px
    style H fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style M fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style N fill:#ffcdd2,stroke:#c62828,stroke-width:2px
EOF

# Module Structure Diagram
cat > "$TEMP_DIR/module_structure.mmd" << 'EOF'
graph TD
    A["🏛️ TSCSMethods.jl<br/>Main Module"] --> B["🔷 core/<br/>2 files"]
    A --> C["🔍 matching/<br/>9 files"]
    A --> D["⚖️ balancing/<br/>7 files"]
    A --> E["📊 estimation/<br/>9 files"]
    A --> F["🚀 advanced/<br/>3 files"]
    A --> G["🛠️ utilities/<br/>7 files"]
    
    B --> B1["📋 types.jl<br/>37 type definitions"]
    B --> B2["🏗️ construction.jl<br/>makemodel() function"]
    
    C --> C1["⚙️ matching_setup.jl"]
    C --> C2["📏 distancing.jl<br/>Mahalanobis distances"]
    C --> C3["🎯 match.jl<br/>Core matching algorithm"]
    C --> C4["📊 ranking.jl<br/>Match ranking system"]
    C --> C5["🔄 retrieve_matches*.jl<br/>3 specialized files"]
    
    D --> D1["⚖️ balancing.jl<br/>Manual balancing"]
    D --> D2["🤖 autobalancing.jl<br/>P-value optimization"]
    D --> D3["📈 meanbalancing.jl<br/>Mean-based strategy"]
    D --> D4["🔄 fullbalancing.jl<br/>Multiple strategies"]
    
    E --> E1["📊 estimation.jl<br/>ATT calculation"]
    E --> E2["🎲 bootstrapping.jl<br/>Statistical inference"]
    E --> E3["📈 overall.jl<br/>Result compilation"]
    E --> E4["👥 estimation_stratified.jl<br/>Group analysis"]
    
    F --> F1["👥 stratification.jl<br/>Subgroup methods"]
    F --> F2["🎯 refine.jl<br/>Match refinement"]
    F --> F3["🔢 groupindices.jl<br/>Index utilities"]
    
    G --> G1["🔍 inspection.jl<br/>Diagnostic tools"]
    G --> G2["💾 storage.jl<br/>Save/load functions"]
    G --> G3["📋 example_data.jl<br/>Sample datasets"]
    
    style A fill:#1976d2,color:#fff,stroke:#0d47a1,stroke-width:3px
    style B fill:#4fc3f7,stroke:#0288d1,stroke-width:2px
    style C fill:#81c784,stroke:#388e3c,stroke-width:2px
    style D fill:#ffb74d,stroke:#f57c00,stroke-width:2px
    style E fill:#f06292,stroke:#c2185b,stroke-width:2px
    style F fill:#ba68c8,stroke:#7b1fa2,stroke-width:2px
    style G fill:#90a4ae,stroke:#546e7a,stroke-width:2px
EOF

# Generate SVG files using mmdc
echo "Converting advanced_workflow.mmd to SVG..."
mmdc -i "$TEMP_DIR/advanced_workflow.mmd" -o "$ASSETS_DIR/advanced_workflow.svg"

echo "Converting data_validation.mmd to SVG..."
mmdc -i "$TEMP_DIR/data_validation.mmd" -o "$ASSETS_DIR/data_validation.svg"

echo "Converting module_structure.mmd to SVG..."
mmdc -i "$TEMP_DIR/module_structure.mmd" -o "$ASSETS_DIR/module_structure.svg"

echo "Additional diagrams generated successfully!"
echo "SVG files saved to: $ASSETS_DIR"