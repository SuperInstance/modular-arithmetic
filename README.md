# modular-arithmetic

A pure-Rust **modular arithmetic toolkit** providing fast exponentiation, the
extended Euclidean algorithm, Chinese Remainder Theorem solver, Euler's totient
function, primitive roots, and discrete logarithms.

## Why This Matters

Modular arithmetic is the foundation of number theory and modern cryptography.
RSA, Diffie-Hellman, DSA, and many post-quantum schemes all rely on operations
in ℤ/nℤ. The Chinese Remainder Theorem enables efficient multi-precision arithmetic.
Euler's totient function connects to the structure of multiplicative groups.
Primitive roots and discrete logarithms underpin Diffie-Hellman key exchange and
the Digital Signature Algorithm.

This library brings all these tools together in one clean, well-documented,
and thoroughly tested package.

## Features

- **Modular exponentiation** — O(log n) via repeated squaring
- **Extended Euclidean Algorithm** — GCD + Bézout coefficients
- **Modular inverse** — solve a·x ≡ 1 (mod m)
- **Chinese Remainder Theorem** — system of congruences → single solution
- **Euler's totient function** φ(n) — count of integers coprime to n
- **Primitive roots** — find the smallest primitive root modulo n
- **Discrete logarithm** — baby-step giant-step algorithm
- **Legendre & Jacobi symbols** — quadratic residuosity

## Mathematical Background

### Modular Exponentiation

Computing a^b mod m efficiently uses **repeated squaring**:

```
a^b = a^(b₀ + 2b₁ + 4b₂ + ...) = a^b₀ · (a²)^b₁ · (a⁴)^b₂ · ...
```

This reduces the number of multiplications from O(b) to O(log b).

### Extended Euclidean Algorithm

For integers a and b, the algorithm finds integers x and y such that:

```
ax + by = gcd(a, b)
```

This directly gives the modular inverse: if gcd(a, m) = 1, then x ≡ a⁻¹ (mod m).

### Chinese Remainder Theorem (CRT)

Given pairwise coprime moduli m₁, m₂, ..., mₖ, the system:

```
x ≡ r₁ (mod m₁)
x ≡ r₂ (mod m₂)
...
x ≡ rₖ (mod mₖ)
```

has a unique solution modulo M = m₁ · m₂ · ... · mₖ.

The CRT is used in RSA to speed up decryption by a factor of ~4 (using CRT-based
private key operations with p and q separately).

### Euler's Totient Function

φ(n) counts integers in [1, n] coprime to n. For prime p: φ(p) = p - 1.
For n = p₁^a₁ · p₂^a₂ · ... : φ(n) = n · ∏(1 - 1/pᵢ).

**Euler's theorem**: If gcd(a, n) = 1, then a^φ(n) ≡ 1 (mod n).
This generalizes Fermat's little theorem (a^(p-1) ≡ 1 mod p for prime p).

### Primitive Roots

A primitive root g modulo n is a generator of the multiplicative group (ℤ/nℤ)*.
Every element of (ℤ/nℤ)* can be written as g^k for some k.

Primitive roots exist iff n = 1, 2, 4, p^k, or 2p^k for odd prime p.

### Discrete Logarithm

Given g (a generator) and h, find x such that g^x ≡ h (mod n).
The **baby-step giant-step** algorithm solves this in O(√n) time and space.

## Usage

```toml
[dependencies]
modular-arithmetic = "0.1.0"
```

### Fast Exponentiation

```rust
use modular_arithmetic::mod_pow;

// 2^100 mod 10^9+7
let result = mod_pow(2, 100, 1_000_000_007);
println!("2^100 mod 10^9+7 = {}", result);
```

### Extended GCD and Inverse

```rust
use modular_arithmetic::{extended_gcd, mod_inverse};

// Bézout coefficients
let (g, x, y) = extended_gcd(35, 15);
assert_eq!(g, 5);
assert_eq!(35 * x + 15 * y, 5);

// Modular inverse
let inv = mod_inverse(3, 7).unwrap(); // 5, since 3*5 = 15 ≡ 1 (mod 7)
```

### Chinese Remainder Theorem

```rust
use modular_arithmetic::chinese_remainder;

// x ≡ 2 (mod 3)
// x ≡ 3 (mod 5)
// x ≡ 2 (mod 7)
// Solution: x = 23 (mod 105)
let (x, m) = chinese_remainder(&[(2, 3), (3, 5), (2, 7)]).unwrap();
assert_eq!(x, 23);
assert_eq!(m, 105);
```

### Euler's Totient and Primitive Roots

```rust
use modular_arithmetic::{euler_totient, primitive_root, is_primitive_root};

// φ(12) = 4 (numbers coprime to 12: {1, 5, 7, 11})
assert_eq!(euler_totient(12), 4);

// Smallest primitive root mod 7
let g = primitive_root(7).unwrap(); // 3
assert!(is_primitive_root(g, 7));
```

### Discrete Logarithm

```rust
use modular_arithmetic::discrete_log;

// Find x such that 3^x ≡ 2 (mod 7)
let x = discrete_log(3, 2, 7);
assert_eq!(x, Some(2));
```

### Legendre and Jacobi Symbols

```rust
use modular_arithmetic::{legendre_symbol, jacobi_symbol};

// Is 2 a quadratic residue mod 7?
assert_eq!(legendre_symbol(2, 7), 1);  // Yes

// Is 3 a quadratic residue mod 7?
assert_eq!(legendre_symbol(3, 7), -1); // No
```

## API Reference

| Function | Description |
|---|---|
| `mod_pow(base, exp, m)` | Fast modular exponentiation |
| `extended_gcd(a, b)` | GCD + Bézout coefficients |
| `mod_inverse(a, m)` | Modular multiplicative inverse |
| `chinese_remainder(pairs)` | CRT solver |
| `euler_totient(n)` | Euler's φ function |
| `is_primitive_root(g, n)` | Check if g is a primitive root |
| `primitive_root(n)` | Find smallest primitive root |
| `discrete_log(g, h, n)` | Baby-step giant-step |
| `legendre_symbol(a, p)` | Legendre symbol (a/p) |
| `jacobi_symbol(a, n)` | Jacobi symbol (a/n) |

## Testing

```bash
cargo test
```

## License

MIT License. See [LICENSE](LICENSE) for details.
