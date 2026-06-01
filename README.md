# secp256k1 Batch Inversion Benchmark

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

High-performance batch inversion optimization for secp256k1 elliptic curve operations, achieving **5-7x speedup** in sequential private key search operations using Montgomery's trick.

## 📊 Key Results

| Method | Speed | Inversions | Speedup |
|--------|-------|------------|---------|
| Naive (Affine) | 2,745 keys/s | 2,000 | 1.00x |
| Projective (Individual) | 2,840 keys/s | 2,000 | 1.03x |
| **Batch Inversion (1000)** | **14,863 keys/s** | **2** | **5.41x** ✅ |
| **Batch Inversion (100)** | **20,542 keys/s** | **50** | **6.88x** ✅ |

**Tested with 2,000 keys on Intel Core i7, Python 3.10**

---

## 🎯 The Discovery

### The Problem

When performing sequential private key searches on the secp256k1 elliptic curve (such as in Bitcoin Puzzle challenges), the primary bottleneck is **field inversion** during coordinate conversion.

Each conversion from projective to affine coordinates requires:

$z⁻¹ = z^(p-2) mod p$

This operation costs approximately **100 field multiplications** using Fermat's little theorem, making it extremely expensive when performed millions of times.

### The Solution

By applying **Montgomery's trick** (batch inversion), we can convert N individual inversions into:
- **1 inversion** of the product
- **~3N multiplications** to extract all inverses

This reduces the cost from `N × 100 multiplications` to `1 × 100 + 3N multiplications`, achieving dramatic speedups for N > 10.

### The Insight

Projective coordinates alone do **not** provide speedup because the final conversion to affine coordinates still requires N inversions. The real optimization comes from:

1. **Computing in projective coordinates** (fast, no inversions)
2. **Batching coordinate conversions** (1 inversion per batch)
3. **Using incremental addition** (Q_{i+1} = Q_i + G instead of computing kG from scratch)

---

## 🧮 Mathematical Foundation

### Montgomery's Trick Algorithm

**Input:** Array of field elements `[z₁, z₂, ..., zₙ]`

**Output:** Array of inverses `[z₁⁻¹, z₂⁻¹, ..., zₙ⁻¹]`

#### Step 1: Prefix Products (N-1 multiplications)

$prefix[0] = z₁$

$prefix[1] = z₁ · z₂$

$prefix[2] = z₁ · z₂ · z₃$

$...$

$prefix[n-1] = z₁ · z₂ · ... · zₙ$


#### Step 2: Single Inversion (~100 multiplications)

$inv = (z₁ · z₂ · ... · zₙ)⁻¹$

$= z₁⁻¹ · z₂⁻¹ · ... · zₙ⁻¹$


Using Fermat's little theorem: `inv = (prefix[n-1])^(p-2) mod p`

#### Step 3: Backward Pass (2N multiplications)

For $i$ from $n-1$ down to $$1:

$result[i] = inv · prefix[i-1]$

This equals: $(z₁⁻¹·...·zₙ⁻¹) · (z₁·...·z_{i-1}) = zᵢ⁻¹$ ✓

$inv = inv · zᵢ$

This equals: $z₁⁻¹·...·z_{i-1}⁻¹$

$result[0] = inv This equals z₁⁻¹$ ✓

### Example: 4 Elements

**Input:** `[z₁, z₂, z₃, z₄]`

**Prefix products:**

prefix = $[z₁, z₁z₂, z₁z₂z₃, z₁z₂z₃z₄]$

Cost: 3 multiplications

**Single inversion:**

$inv = (z₁z₂z₃z₄)⁻¹$

Cost: ~100 multiplications

**Backward pass:**

$result[3] = inv · prefix[2]$

$= (z₁z₂z₃z₄)⁻¹ · (z₁z₂z₃)$

$= z₄⁻¹ $ ✓

$inv = inv · z₄ = z₁⁻¹z₂⁻¹z₃⁻¹$

$result[2] = inv · prefix[1]$

$= (z₁⁻¹z₂⁻¹z₃⁻¹) · (z₁z₂)$

$= z₃⁻¹$ ✓

$inv = inv · z₃ = z₁⁻¹z₂⁻¹$ 

$result[1] = inv · prefix[0]$

$= (z₁⁻¹z₂⁻¹) · z₁$

$= z₂⁻¹ $ ✓

$inv = inv · z₂ = z₁⁻¹$

$result[0] = inv = z₁⁻¹ $✓

Cost: 8 multiplications


**Total cost:** 3 + 100 + 8 = **111 multiplications**

**Naive approach:** 4 × 100 = **400 multiplications**

**Speedup:** 400 / 111 = **3.6x**

### Scaling Analysis

| N (elements) | Naive Cost | Montgomery Cost | Speedup |
|--------------|------------|-----------------|---------|
| 10 | 1,000 | 100 + 30 + 20 = 150 | 6.7x |
| 100 | 10,000 | 100 + 300 + 200 = 600 | 16.7x |
| 1,000 | 100,000 | 100 + 3,000 + 2,000 = 5,100 | 19.6x |
| 10,000 | 1,000,000 | 100 + 30,000 + 20,000 = 50,100 | 20.0x |

**Asymptotic speedup:** ~20x for large N

---

## 🚀 Quick Start

### Installation

```bash
git clone https://Anvexis/secp256k1-batch-inversion.git
cd secp256k1-batch-inversion
python3 test_batch_inversion.py
```

Expected Output
```python
======================================================================
SECP256K1 BATCH INVERSION BENCHMARK
Testing: Montgomery's Trick + Incremental Addition
======================================================================

Number of keys to test: 2,000
This may take a few minutes...

======================================================================
Benchmark 1: Naive Approach (Affine Coordinates)
======================================================================
Time: 0.729 seconds
Speed: 2,745 keys/s
Inversions: 2,000

======================================================================
Benchmark 2: Projective Coordinates (Individual Conversion)
======================================================================
Time: 0.704 seconds
Speed: 2,840 keys/s
Inversions: 2,000

======================================================================
Benchmark 3: Batch Inversion (batch size: 1000)
======================================================================
Time: 0.135 seconds
Speed: 14,863 keys/s
Inversions: 2

======================================================================
METHOD COMPARISON
======================================================================
Method                                   Speed       Speedup
----------------------------------------------------------------------
Naive (Affine)                            2,745 keys/s    1.00x
Projective (Individual)                   2,840 keys/s    1.03x
Batch Inversion (1000)                   14,863 keys/s    5.41x

======================================================================
CONCLUSIONS
======================================================================
✅ Batch Inversion provides 5.23x speedup compared to
   individual projective coordinate conversion
✅ Optimal batch size: 100-1000 points
✅ Montgomery's trick reduces inversions: 2000 → 2
✅ Real-world impact: Bitcoin Puzzle #71 search time reduced from
   13,300 years to 2,500 years (single RTX 4090 GPU)
======================================================================

```

### C++ Implementation (Production)
### For production use, implement in C++ with optimizations:

```cpp
std::vector<uint256_t> batch_inversion(const std::vector<uint256_t>& elements) {
    size_t n = elements.size();
    std::vector<uint256_t> prefix(n), result(n);
    
    // Prefix products
    prefix[0] = elements[0];
    for (size_t i = 1; i < n; ++i) {
        prefix[i] = mod_mul(prefix[i-1], elements[i]);
    }
    
    // Single inversion
    uint256_t inv = mod_inverse(prefix[n-1]);
    
    // Backward pass
    for (size_t i = n - 1; i > 0; --i) {
        result[i] = mod_mul(inv, prefix[i-1]);
        inv = mod_mul(inv, elements[i]);
    }
    result[0] = inv;
    
    return result;
}
```

### CUDA Implementation (GPU)
### For GPU acceleration:

```cpp
__global__ void batch_inversion_kernel(
    const uint256_t* z_values,
    uint256_t* z_invs,
    size_t batch_size
) {
    extern __shared__ uint256_t shared_prefix[];
    
    size_t tid = threadIdx.x;
    
    // Load into shared memory
    shared_prefix[tid] = z_values[blockIdx.x * batch_size + tid];
    __syncthreads();
    
    // Parallel prefix sum (work-efficient algorithm)
    // ... (implementation details)
    
    // Single inversion by thread 0
    __shared__ uint256_t total_inv;
    if (tid == 0) {
        total_inv = mod_inverse(shared_prefix[batch_size - 1]);
    }
    __syncthreads();
    
    // Backward pass
    // ... (implementation details)
}
```

### Multi-GPU Distribution

```
def distributed_search(gpu_count, batch_size=1000):
    total_keys = 2**70
    keys_per_gpu = total_keys // gpu_count
    
    processes = []
    for gpu_id in range(gpu_count):
        start_key = 2**70 + gpu_id * keys_per_gpu
        p = Process(target=search_on_gpu, 
                   args=(gpu_id, start_key, keys_per_gpu, batch_size))
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()
```

## What This Project Does NOT Claim

This project presents an engineering optimization for elliptic curve computations on secp256k1.

The batch inversion technique (Montgomery's Trick):

* Reduces the number of field inversions required during coordinate conversion.
* Improves throughput of sequential elliptic curve point generation.
* Accelerates ECC workloads used in cryptographic research and benchmarking.

This project does **not**:

* Break secp256k1.
* Solve the Elliptic Curve Discrete Logarithm Problem (ECDLP).
* Reduce the asymptotic complexity of Pollard Rho, Pollard Kangaroo, or Baby-Step Giant-Step.
* Provide a practical attack against Bitcoin private keys.

The achieved speedup is a constant-factor engineering improvement rather than a cryptanalytic breakthrough.

## Practical Applications

The techniques demonstrated in this repository are applicable to:

* Pollard Rho implementations
* Pollard Kangaroo implementations
* Baby-Step Giant-Step (BSGS) solvers
* Vanity address generation
* secp256k1 research frameworks
* GPU ECC computation pipelines
* Large-scale point table generation
* Cryptographic benchmarking

Batch inversion is particularly beneficial when thousands of projective points must be converted to affine coordinates simultaneously.

## Limitations

Benchmark results presented in this repository were obtained using Python implementations and should not be interpreted as theoretical speedups.

Observed acceleration depends on:

* CPU architecture
* Compiler optimizations
* Field arithmetic implementation
* Batch size
* Memory bandwidth
* Cache efficiency

For production-grade C/C++ implementations using optimized secp256k1 arithmetic, absolute performance figures may differ significantly.

The theoretical benefit of Montgomery's trick approaches approximately 20× for very large batches, but real-world implementations typically achieve smaller gains due to surrounding computation costs.

## Reproducibility

The benchmark results presented in this repository can be reproduced using the included test suite.

### Environment

Recommended environment:

* Python 3.8+
* Linux (Ubuntu 22.04 tested)
* GCC 11+
* Intel or AMD x86_64 CPU
* 16 GB RAM or more

### Clone Repository

```bash
git clone https://github.com/Anvexis/secp256k1-batch-inversion.git
cd secp256k1-batch-inversion
```

### Run Benchmark

```bash
python3 test_batch_inversion.py
```

### Expected Output

The benchmark will evaluate:

1. Affine coordinate computation
2. Jacobian coordinates with individual conversion
3. Jacobian coordinates with Montgomery batch inversion
4. Jacobian coordinates with batch inversion and incremental addition

Example output:

```text
Naive (Affine):                      ~2,700 keys/s
Jacobian (Individual Conversion):    ~2,800 keys/s
Jacobian + Batch Inversion:         ~15,000 keys/s
Jacobian + Batch + Incremental:     ~20,000 keys/s
```

Exact results depend on:

* CPU architecture
* Compiler version
* Python version
* Cache size
* Batch size

### Verification

To verify correctness:

* All generated public keys must match the affine reference implementation.
* Batch inversion outputs must match individual modular inversions.
* Hash160 results must remain identical across all methods.

The optimization changes only performance characteristics and does not alter cryptographic correctness.

## Performance Comparison

The following table summarizes the tested approaches.

| Method                                            | Point Generation | Coordinate System | Inversions  | Relative Speed |
| ------------------------------------------------- | ---------------- | ----------------- | ----------- | -------------- |
| Affine                                            | kG from scratch  | Affine            | N           | 1.0x           |
| Jacobian                                          | kG from scratch  | Jacobian          | N           | 1.0–1.2x       |
| Jacobian + Batch Inversion                        | kG from scratch  | Jacobian          | 1 per batch | 5–7x           |
| Jacobian + Batch Inversion + Incremental Addition | Q(i+1)=Q(i)+G    | Jacobian          | 1 per batch | 7–20x*         |

* Depends heavily on implementation details, batch size, hardware, and surrounding cryptographic workload.

### Interpretation

#### Affine Coordinates

Each point operation requires modular inversion.

Advantages:

* Simple implementation
* Easy debugging

Disadvantages:

* Inversion-heavy
* Poor scalability

#### Jacobian Coordinates

Field inversions are delayed until final conversion.

Advantages:

* Faster point addition
* Standard ECC optimization

Disadvantages:

* Final affine conversion still requires N inversions

#### Jacobian + Batch Inversion

Uses Montgomery's trick to replace N inversions with a single inversion.

Advantages:

* Significant reduction in inversion cost
* Excellent for large batches

Disadvantages:

* Requires temporary storage
* Benefits appear only when processing many points

#### Jacobian + Batch Inversion + Incremental Addition

Combines:

1. Jacobian arithmetic
2. Montgomery batch inversion
3. Sequential point generation

Instead of computing:

kG

for every key independently, the next point is generated as:

Q(i+1) = Q(i) + G

This avoids repeated scalar multiplication and is the preferred strategy for large-scale sequential key scanning.

Advantages:

* Maximum throughput
* Minimal inversion count
* Excellent GPU compatibility

Disadvantages:

* Applicable primarily to sequential key ranges
* Requires careful state management
