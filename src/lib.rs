//! # Modular Arithmetic Toolkit
//!
//! A pure-Rust library for modular arithmetic operations, including fast
//! modular exponentiation, the extended Euclidean algorithm, Chinese Remainder
//! Theorem, Euler's totient function, primitive roots, and discrete logarithms.

/// Modular exponentiation via repeated squaring.
///
/// Computes base^exp mod modulus in O(log exp) time.
///
/// # Panics
/// Panics if modulus ≤ 0.
pub fn mod_pow(base: i64, exp: i64, modulus: i64) -> i64 {
    assert!(modulus > 0, "Modulus must be positive");
    if modulus == 1 { return 0; }
    let mut result = 1i64;
    let mut base = base.rem_euclid(modulus);
    let mut exp = exp;
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % modulus;
        }
        exp >>= 1;
        if exp > 0 {
            base = (base * base) % modulus;
        }
    }
    result
}

/// Extended Euclidean Algorithm.
///
/// Returns (gcd, x, y) such that a*x + b*y = gcd(a, b).
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        return (b, 0, 1);
    }
    let (g, x, y) = extended_gcd(b % a, a);
    (g, y - (b / a) * x, x)
}

/// Modular multiplicative inverse.
///
/// Returns x such that a*x ≡ 1 (mod m), or None if gcd(a, m) ≠ 1.
pub fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(a.rem_euclid(m), m);
    if g != 1 {
        return None;
    }
    Some(x.rem_euclid(m))
}

/// Chinese Remainder Theorem.
///
/// Given pairs (r₁, m₁), (r₂, m₂), ..., finds x such that:
/// - x ≡ r₁ (mod m₁)
/// - x ≡ r₂ (mod m₂)
/// - ...
///
/// Returns (x, M) where M = m₁ * m₂ * ... and x is the unique solution mod M.
/// Returns None if the moduli are not pairwise coprime.
pub fn chinese_remainder(congruences: &[(i64, i64)]) -> Option<(i64, i64)> {
    if congruences.is_empty() {
        return None;
    }

    let mut x = congruences[0].0.rem_euclid(congruences[0].1);
    let mut m = congruences[0].1;

    for &(ri, mi) in &congruences[1..] {
        let (g, u, _) = extended_gcd(m, mi);
        if (ri - x) % g != 0 {
            return None; // No solution
        }
        let lcm = m / g * mi;
        let diff = ri - x;
        let inv = (diff / g).rem_euclid(mi);
        x = (x + m * ((u * inv) % (mi / g))).rem_euclid(lcm);
        if x < 0 { x += lcm; }
        m = lcm;
    }

    Some((x, m))
}

/// Euler's totient function φ(n).
///
/// Counts the number of integers 1 ≤ k ≤ n that are coprime to n.
pub fn euler_totient(n: i64) -> i64 {
    if n <= 0 { return 0; }
    if n == 1 { return 1; }

    let mut result = n;
    let mut m = n;

    let mut p = 2;
    while p * p <= m {
        if m % p == 0 {
            while m % p == 0 {
                m /= p;
            }
            result -= result / p;
        }
        p += 1;
    }
    if m > 1 {
        result -= result / m;
    }

    result
}

/// Check if g is a primitive root modulo n.
///
/// A primitive root generates the entire multiplicative group modulo n.
/// g is a primitive root iff for every prime factor q of φ(n): g^(φ(n)/q) ≢ 1 (mod n).
pub fn is_primitive_root(g: i64, n: i64) -> bool {
    if n <= 2 { return g % n == 1; }

    let g = g.rem_euclid(n);
    if g <= 1 { return false; }

    // Check gcd(g, n) = 1
    let (gcd_val, _, _) = extended_gcd(g, n);
    if gcd_val != 1 { return false; }

    let phi = euler_totient(n);

    // Get prime factors of φ(n)
    let factors = prime_factors_unique(phi);

    // Check g^(φ/q) ≢ 1 for each prime factor q
    for q in &factors {
        if mod_pow(g, phi / q, n) == 1 {
            return false;
        }
    }

    true
}

/// Find the smallest primitive root modulo n.
///
/// Returns None if no primitive root exists (i.e., n is not 1, 2, 4, p^k, or 2p^k
/// for odd prime p).
pub fn primitive_root(n: i64) -> Option<i64> {
    if n <= 1 { return None; }
    if n == 2 { return Some(1); }
    if n == 4 { return Some(3); }

    // Check if primitive roots exist: n = p^k or 2p^k for odd prime p
    let factors = prime_factors_with_multiplicity(n);
    if factors.len() == 1 && factors[0].0 != 2 {
        // n = p^k, OK
    } else if factors.len() == 2 && factors[0].0 == 2 && factors[0].1 == 1 {
        // n = 2 * p^k, OK
    } else if factors.len() == 1 && factors[0].0 == 2 && factors[0].1 >= 3 {
        return None; // 2^k for k≥3 has no primitive roots
    } else if factors.len() > 2 {
        return None;
    }

    let phi = euler_totient(n);
    let prime_factors = prime_factors_unique(phi);

    for g in 2..n {
        let mut is_root = true;
        let (gcd_val, _, _) = extended_gcd(g, n);
        if gcd_val != 1 { continue; }

        for q in &prime_factors {
            if mod_pow(g, phi / q, n) == 1 {
                is_root = false;
                break;
            }
        }

        if is_root {
            return Some(g);
        }
    }

    None
}

/// Discrete logarithm via baby-step giant-step.
///
/// Solves g^x ≡ h (mod n) for x. Returns None if no solution.
pub fn discrete_log(g: i64, h: i64, n: i64) -> Option<i64> {
    let g = g.rem_euclid(n);
    let h = h.rem_euclid(n);

    if g == 0 && h == 0 { return Some(1); }
    if g == 0 { return None; }
    if h == 1 { return Some(0); }
    if g == h { return Some(1); }

    let n_sqrt = ((n as f64).sqrt() as i64) + 1;

    // Baby steps: store g^j → j
    let mut table = std::collections::HashMap::new();
    let mut power = 1i64;
    for j in 0..=n_sqrt {
        table.entry(power).or_insert(j);
        power = (power * g) % n;
    }

    // Giant step: g^(-m)
    let g_inv = mod_inverse(g, n)?;
    let g_m_inv = mod_pow(g_inv, n_sqrt, n);

    let mut gamma = h;
    for i in 0..=n_sqrt {
        if let Some(&j) = table.get(&gamma) {
            let x = (i * n_sqrt + j) % (euler_totient(n).max(1));
            if mod_pow(g, x, n) == h {
                return Some(x);
            }
            // Try the raw value
            let raw_x = i * n_sqrt + j;
            if mod_pow(g, raw_x, n) == h {
                return Some(raw_x);
            }
        }
        gamma = (gamma * g_m_inv) % n;
    }

    None
}

/// Get unique prime factors of n.
fn prime_factors_unique(mut n: i64) -> Vec<i64> {
    let mut factors = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            factors.push(d);
            while n % d == 0 {
                n /= d;
            }
        }
        d += 1;
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

/// Get prime factors with multiplicity: (p, k) where p^k | n.
fn prime_factors_with_multiplicity(mut n: i64) -> Vec<(i64, i64)> {
    let mut factors = Vec::new();
    let mut d = 2;
    while d * d <= n {
        if n % d == 0 {
            let mut count = 0;
            while n % d == 0 {
                n /= d;
                count += 1;
            }
            factors.push((d, count));
        }
        d += 1;
    }
    if n > 1 {
        factors.push((n, 1));
    }
    factors
}

/// Compute the Legendre symbol (a/p) for odd prime p.
///
/// Returns 1 if a is a QR, -1 if a is a QNR, 0 if p divides a.
pub fn legendre_symbol(a: i64, p: i64) -> i64 {
    let a = a.rem_euclid(p);
    if a == 0 { return 0; }
    let ls = mod_pow(a, (p - 1) / 2, p);
    if ls == p - 1 { -1 } else { ls }
}

/// Compute the Jacobi symbol (a/n) for odd positive n.
pub fn jacobi_symbol(mut a: i64, mut n: i64) -> i64 {
    assert!(n > 0 && n % 2 == 1, "n must be a positive odd integer");
    a = a.rem_euclid(n);
    let mut result = 1;

    while a != 0 {
        while a % 2 == 0 {
            a /= 2;
            if n % 8 == 3 || n % 8 == 5 {
                result = -result;
            }
        }
        std::mem::swap(&mut a, &mut n);
        if a % 4 == 3 && n % 4 == 3 {
            result = -result;
        }
        a %= n;
    }

    if n == 1 { result } else { 0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mod_pow_basic() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 0, 7), 1);
        assert_eq!(mod_pow(5, 1, 7), 5);
    }

    #[test]
    fn test_mod_pow_large() {
        assert_eq!(mod_pow(2, 100, 1000000007), 976371285);
    }

    #[test]
    fn test_extended_gcd() {
        let (g, x, y) = extended_gcd(35, 15);
        assert_eq!(g, 5);
        assert_eq!(35 * x + 15 * y, 5);
    }

    #[test]
    fn test_extended_gcd_coprime() {
        let (g, x, y) = extended_gcd(17, 12);
        assert_eq!(g, 1);
        assert_eq!(17 * x + 12 * y, 1);
    }

    #[test]
    fn test_mod_inverse_exists() {
        let inv = mod_inverse(3, 7).unwrap();
        assert_eq!((3 * inv) % 7, 1);
        assert_eq!(inv, 5);
    }

    #[test]
    fn test_mod_inverse_not_exists() {
        assert!(mod_inverse(2, 4).is_none());
    }

    #[test]
    fn test_chinese_remainder_simple() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
        let result = chinese_remainder(&[(2, 3), (3, 5), (2, 7)]);
        let (x, m) = result.unwrap();
        assert_eq!(m, 105);
        assert_eq!(x % 3, 2);
        assert_eq!(x % 5, 3);
        assert_eq!(x % 7, 2);
    }

    #[test]
    fn test_chinese_remainder_two() {
        // x ≡ 1 (mod 2), x ≡ 2 (mod 3)
        let result = chinese_remainder(&[(1, 2), (2, 3)]);
        let (x, m) = result.unwrap();
        assert_eq!(m, 6);
        assert_eq!(x, 5);
    }

    #[test]
    fn test_euler_totient_prime() {
        assert_eq!(euler_totient(7), 6);
        assert_eq!(euler_totient(13), 12);
    }

    #[test]
    fn test_euler_totient_composite() {
        assert_eq!(euler_totient(12), 4); // {1, 5, 7, 11}
        assert_eq!(euler_totient(10), 4); // {1, 3, 7, 9}
        assert_eq!(euler_totient(1), 1);
    }

    #[test]
    fn test_euler_totient_prime_power() {
        assert_eq!(euler_totient(9), 6);  // φ(3²) = 3² - 3 = 6
        assert_eq!(euler_totient(8), 4);  // φ(2³) = 2³ - 2² = 4
    }

    #[test]
    fn test_is_primitive_root() {
        // 3 is a primitive root mod 7
        assert!(is_primitive_root(3, 7));
        // 2 is NOT a primitive root mod 7
        assert!(!is_primitive_root(2, 7));
    }

    #[test]
    fn test_primitive_root_finds_smallest() {
        // Smallest primitive root mod 7 is 3
        assert_eq!(primitive_root(7), Some(3));
        // Smallest primitive root mod 11 is 2
        assert_eq!(primitive_root(11), Some(2));
    }

    #[test]
    fn test_discrete_log_basic() {
        // 3^x ≡ 2 mod 7 → x = 2
        assert_eq!(discrete_log(3, 2, 7), Some(2));
    }

    #[test]
    fn test_discrete_log_identity() {
        assert_eq!(discrete_log(3, 1, 7), Some(0));
    }

    #[test]
    fn test_legendre_symbol() {
        assert_eq!(legendre_symbol(2, 7), 1);   // 2 is QR mod 7
        assert_eq!(legendre_symbol(3, 7), -1);   // 3 is QNR mod 7
        assert_eq!(legendre_symbol(7, 7), 0);     // 7 ≡ 0 mod 7
    }

    #[test]
    fn test_jacobi_symbol() {
        assert_eq!(jacobi_symbol(2, 7), 1);
        assert_eq!(jacobi_symbol(5, 3), -1);
    }

    #[test]
    fn test_fermat_little_theorem() {
        // a^(p-1) ≡ 1 (mod p) for prime p not dividing a
        for a in [2, 3, 5] {
            assert_eq!(mod_pow(a, 12, 13), 1);
        }
    }

    #[test]
    fn test_euler_theorem() {
        // a^φ(n) ≡ 1 (mod n) for gcd(a,n) = 1
        assert_eq!(mod_pow(7, euler_totient(15), 15), 1);
    }
}
