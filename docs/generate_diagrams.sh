#!/bin/bash

# Script to generate SVG diagrams from Mermaid source files
# Requires: npm install -g @mermaid-js/mermaid-cli

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSETS_DIR="$SCRIPT_DIR/src/assets/images"

# Create assets directory if it doesn't exist
mkdir -p "$ASSETS_DIR"

# Create temporary directory for mermaid files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Generating diagrams..."

# Methodology Flow Diagram
cat > "$TEMP_DIR/methodology_flow.mmd" << 'EOF'
flowchart TD
    A[Time-Series Cross-Sectional Data] --> B[Treatment Event Detection]
    B --> C[Define Time Windows]
    C --> D["F: Post-Treatment Periods<br/>(e.g., 1:10)"]
    C --> E["L: Pre-Treatment Periods<br/>(e.g., -20:-1)"]
    
    D --> F[Matching Process]
    E --> F
    
    F --> G[Mahalanobis Distance Calculation]
    G --> H[Covariate Balancing]
    H --> I[Treatment Effect Estimation]
    I --> J[Bootstrap Inference]
    J --> K["ATT with Confidence Intervals<br/>📊 Final Results"]
    
    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style K fill:#c8e6c9,stroke:#2e7d32,stroke-width:3px
    style F fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style H fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
EOF

# User Workflow Diagram
cat > "$TEMP_DIR/user_workflow.mmd" << 'EOF'
flowchart TD
    A["📥 Load Data<br/>example_data()"] --> B["🏗️ Create Model<br/>makemodel()"]
    B --> C["⏰ Specify Time Windows<br/>F: post-treatment<br/>L: pre-treatment"]
    C --> D["📋 Define Covariates<br/>& Time-Varying Dict"]
    
    D --> E["🔍 Matching Phase"]
    E --> F["match!(model, data)"]
    F --> G["⚖️ Balancing Phase"]
    G --> H["balance!(model, data)"]
    H --> I["📊 Estimation Phase"]
    I --> J["estimate!(model, data)"]
    
    J --> K["📈 Results Available"]
    K --> L["model.overall.ATT"]
    K --> M["Confidence Intervals<br/>p05, p95"]
    K --> N["Bootstrap Distributions"]
    
    style A fill:#e3f2fd,stroke:#1565c0,stroke-width:2px
    style E fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style G fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
    style I fill:#e8f5e8,stroke:#2e7d32,stroke-width:2px
    style K fill:#fff8e1,stroke:#f57f17,stroke-width:3px
EOF

# Validation Framework Diagram
cat > "$TEMP_DIR/validation_framework.mmd" << 'EOF'
flowchart LR
    A[Core Estimation] --> B{Advanced Features}
    
    B --> C[Stratified Estimation]
    B --> D[Caliper Matching]
    B --> E[Refinement Procedures]
    B --> F[Auto-balancing]
    
    C --> G[Statistical Validation Framework]
    D --> G
    E --> G
    F --> G
    
    G --> H["📊 Coverage Tests<br/>Null DGP"]
    G --> I["🎲 Placebo Tests<br/>Permutation"]
    G --> J["🎯 Noiseless Recovery<br/>Exact Tests"]
    
    H --> K["✅ 96% Coverage<br/>(Target: 93-97%)"]
    I --> L["✅ 6.87% Type I Error<br/>(Target: 3-7%)"]
    J --> M["✅ Exact Recovery<br/>(1e-3 tolerance)"]
    
    style G fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style K fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style L fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    style M fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
EOF

# Package Architecture Diagram
cat > "$TEMP_DIR/package_architecture.mmd" << 'EOF'
graph TB
    subgraph "Core System"
        A[types.jl<br/>Model Types & Data Structures]
        B[construction.jl<br/>makemodel()]
        C[dependencies.jl<br/>External Dependencies]
        D[utils.jl<br/>Utility Functions]
    end
    
    subgraph "Matching System"
        E[matching_setup.jl<br/>Match Preparation]
        F[matching_algo.jl<br/>Core Matching Logic]
        G[matching_utils.jl<br/>Distance & Ranking]
        H[getmatches.jl<br/>Match Retrieval]
    end
    
    subgraph "Balancing System"
        I[balancing.jl<br/>Manual Balance Assessment]
        J[balancing_auto.jl<br/>Auto-balancing Logic]
        K[balancing_utils.jl<br/>Balance Calculations]
    end
    
    subgraph "Estimation System"
        L[estimation.jl<br/>Treatment Effect Estimation]
        M[bootstrap.jl<br/>Bootstrap Inference]
        N[bayesfactor.jl<br/>Evidence Quantification]
    end
    
    subgraph "Advanced Features"
        O[stratified.jl<br/>Stratified Analysis]
        P[refinement.jl<br/>Match Refinement]
        Q[calipers.jl<br/>Caliper Matching]
        R[groupindices.jl<br/>Group Management]
    end
    
    subgraph "Utilities"
        S[data_utils.jl<br/>Data Processing]
        T[example_data.jl<br/>Sample Datasets]
        U[simulate_tscs.jl<br/>Data Generation]
        V[validation.jl<br/>Statistical Tests]
    end
    
    A --> B
    B --> E
    E --> F
    F --> I
    I --> L
    L --> M
    
    C --> B
    D --> F
    G --> F
    H --> F
    K --> I
    J --> I
    N --> L
    
    O --> L
    P --> F
    Q --> F
    R --> E
    
    S --> B
    T --> S
    U --> S
    V --> L
    
    style A fill:#e3f2fd,stroke:#1565c0,stroke-width:2px
    style L fill:#e8f5e8,stroke:#2e7d32,stroke-width:2px
    style M fill:#fff3e0,stroke:#ef6c00,stroke-width:2px
    style O fill:#f3e5f5,stroke:#7b1fa2,stroke-width:2px
EOF

# Generate SVG files using mmdc
echo "Converting methodology_flow.mmd to SVG..."
mmdc -i "$TEMP_DIR/methodology_flow.mmd" -o "$ASSETS_DIR/methodology_flow.svg"

echo "Converting user_workflow.mmd to SVG..."
mmdc -i "$TEMP_DIR/user_workflow.mmd" -o "$ASSETS_DIR/user_workflow.svg"

echo "Converting validation_framework.mmd to SVG..."
mmdc -i "$TEMP_DIR/validation_framework.mmd" -o "$ASSETS_DIR/validation_framework.svg"

echo "Converting package_architecture.mmd to SVG..."
mmdc -i "$TEMP_DIR/package_architecture.mmd" -o "$ASSETS_DIR/package_architecture.svg"

echo "All diagrams generated successfully!"
echo "SVG files saved to: $ASSETS_DIR"