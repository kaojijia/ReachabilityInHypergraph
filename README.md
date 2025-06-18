

# HyperReachability Benchmark Suite

This repository contains benchmark implementations for various hypergraph reachability algorithms, including both static and dynamic query methods.

## Overview

The benchmark suite includes three main test programs:

1. **Static Query Benchmark** - Compares four different reachability methods
2. **Dynamic HRTree Benchmark** - Tests dynamic updates on tree-based hypergraph index
3. **Dynamic Weighted Reachability (DWR) Benchmark** - Tests dynamic updates on bottleneck graph index

## Compilation

### Prerequisites
- CMake 3.10 or higher
- C++17 compatible compiler (GCC 7+ or Clang 5+)
- Google Test framework (included as submodule)
- Intel TBB library

### Build Instructions

```bash
# Clone the repository
git clone <repository-url>
cd ReachabilityInHypergraph


# Create build directory
mkdir build
cd build

# Configure and build
cmake ..
make -j$(nproc)
```

The compiled binaries will be located in:
- `build/test/bench_static_query`
- `build/test/bench_dynamic_HIndex`
- `build/test/bench_dynamic_DWR`

## Test Programs

### 1. Static Query Benchmark (`bench_static_query`)

Compares the performance of four hypergraph reachability methods:

#### Usage
```bash
./build/test/bench_static_query --dataset <hypergraph_file>
```

#### Parameters
- `--dataset <filepath>`: Path to the hypergraph dataset file

#### Example
```bash
./build/test/bench_static_query --dataset /path/to/hypergraph.txt
```


### 2. Dynamic HyperIndex Benchmark (`bench_dynamic_HIndex`)

Tests dynamic update performance on hypergraph tree index structure.

#### Usage
```bash
./build/test/bench_dynamic_HIndex --dataset <hypergraph_file>
```

#### Parameters
- `--dataset <filepath>`: Path to the hypergraph dataset file

#### Example
```bash
./build/test/bench_dynamic_HIndex --dataset /path/to/hypergraph.txt
```


### 3. Dynamic Weighted Reachability Benchmark (`bench_dynamic_DWR`)

Tests dynamic update performance on bottleneck graph (DLCR-like) index structure.

#### Usage
```bash
./build/test/bench_dynamic_DWR --dataset <hypergraph_file> 
```

#### Parameters
- `--dataset <filepath>`: Path to the original hypergraph dataset file

#### Example
```bash
./build/test/bench_dynamic_DWR --dataset /path/to/hypergraph
```


## Dataset Format

All benchmarks expect hypergraph files in the following format:
- Each line represents one hyperedge
- Vertices in a hyperedge are separated by spaces
- Vertex IDs should be non-negative integers
- Example:
  ```
  0 1 2
  1 3 4
  2 4 5
  ```

## Output Interpretation

### Static Benchmark Output
- **Build Time**: Time to construct each index type
- **Memory Usage**: RAM consumption of each index
- **Query Time**: Average time per query in microseconds
- Results are grouped by k values (minimum intersection size)

### Dynamic Benchmark Output
- **Update Time**: Average time per dynamic operation
- **Operation Mix Analysis**: Projected performance for workloads with different add/remove ratios

