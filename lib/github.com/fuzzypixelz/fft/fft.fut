-- | An implementation of the parallel-FFT algorithm for evaluating a DFT in parallel.

import "../../diku-dk/complex/complex"
import "../../diku-dk/cpprandom/random"

module complex = mk_complex f64
module dist = uniform_real_distribution f64 minstd_rand

-- | Computes the `d`-th root of unity of the ring of complex numbers.
def root (n: i64) = 
    complex.mk_im (2f64 * f64.pi / f64.i64 n) 
    |> complex.neg 
    |> complex.exp

-- | Computes (naively) the `k`th component of the DFT of the `d`-point array `a`.
def dft [n] (a: [n]complex.complex) (k: i64) =
    -- Enumerate the elements of the array using i64 indices
    let a' = zip (0i64..<n) a
    -- The d-th root of unity to power of `j` and the sum index `i`
    let root' i = root n complex.** complex.mk_re (f64.i64 (i * k))
    in map (\(i, z) -> z complex.* root' i) a' |> complex.sum

def dft' [n] (a: [n]complex.complex) = map (dft a) (0..<n)

type^ factorize = i64 -> (i64, i64)

-- | Computes the DFT of the `a` array using the parallel-FFT method,
-- the `factor` function is used to decompose `n` into two integers: `p_max` and `q_max`
def pdft [n] (f: factorize) (a: [n]complex.complex) =
    let (p_max, q_max) = f n
    -- Calculate in parallel p_max DFTs with q_max points
    let aslices = 
        map (\p -> 
            map (\q -> 
                a[q * p_max + p]) 
            (0..<q_max) |> dft') 
        (0..<p_max)
    -- Calculate in parallel q_max DFTs with p_max points
    let root' p q = root n complex.** (complex.mk_re (f64.i64 (p * q)))
    let zslices = 
        map (\q -> 
            map (\p -> 
                root' p q complex.* aslices[p][q]) 
            (0..<p_max) |> dft') 
        (0..<q_max)
    -- Restructure the final transform 
    in map (\k -> let (p, q) = (k / q_max, k % q_max) in zslices[q][p]) (0..<n)

def factorize_eq (n: i64) = let r = f64.i64 n |> f64.sqrt |> i64.f64 in (r, r)

def mk_array (n: i64) =
    let (_, result) = 
        loop (rng, rs) = (minstd_rand.rng_from_seed [42], []) for _i < n do
            let (rng, re) = dist.rand (1, 10) rng
            let (rng, im) = dist.rand (1, 10) rng
            let z = complex.mk re im
            in (rng, rs ++ [z])
    in result
