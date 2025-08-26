# API Reference

This page provides comprehensive documentation for all TSCSMethods.jl functions and types.

## Package Architecture

The following diagram shows the modular architecture of TSCSMethods.jl:

![Package Architecture](assets/images/package_architecture.svg)

The architecture follows a modular design with clear separation of concerns:

- **Core System**: Fundamental types and model construction
- **Matching System**: Distance-based matching algorithms with sophisticated control unit selection
- **Balancing System**: Covariate balance assessment and optimization
- **Estimation System**: Treatment effect estimation with bootstrap inference
- **Advanced Features**: Specialized methods for complex scenarios
- **Utilities**: Data processing, examples, and validation tools

## Complete Function and Type Reference

```@autodocs
Modules = [TSCSMethods]
```