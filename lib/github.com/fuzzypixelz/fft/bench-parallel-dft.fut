import "fft"

-- ==
-- compiled input { 4 }
-- compiled input { 5 }
-- compiled input { 6 }

def main (n: i32) =
    let size = 2 ** (2 * n) |> i64.i32
    let a = mk_array size
    in pdft factorize_eq a