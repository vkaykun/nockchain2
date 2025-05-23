use alloc::vec::Vec;
use crate::stark::field::{fft_in_place, ifft_in_place, Belt};

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__compute_codeword_commitments(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-codeword-commitments
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__compute_lde(ptr: *mut u64) {
    // Safety: VM guarantees pointer is valid and aligned for u64
    if ptr.is_null() {
        return;
    }
    unsafe {
        let len = *ptr as usize;
        let log_n = *ptr.add(1) as usize;
        let coeffs = core::slice::from_raw_parts_mut(ptr.add(2), len);

        // Copy to Belt vector
        let mut vec: Vec<Belt> = Vec::with_capacity(len);
        vec.extend(coeffs.iter().map(|&c| Belt(c)));

        fft_in_place(&mut vec, log_n);

        // Write back results
        for (dest, v) in coeffs.iter_mut().zip(vec.into_iter()) {
            *dest = v.0;
        }
    }
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__compute_composition_poly(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-composition-poly
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__compute_deep(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-deep
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__precompute_ntts(_ptr: *mut u64) {
    // TODO: implement native accelerator for +precompute-ntts
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__build_compute_queue(_ptr: *mut u64) {
    // TODO: implement native accelerator for +build-compute-queue
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__build_jute_list(_ptr: *mut u64) {
    // TODO: implement native accelerator for +build-jute-list
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__indirect_to_bits(_ptr: *mut u64) {
    // TODO: implement native accelerator for +indirect-to-bits
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__noun_get_zero_mults(_ptr: *mut u64) {
    // TODO: implement native accelerator for +noun-get-zero-mults
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__pbkdf(_ptr: *mut u64) {
    // TODO: implement native accelerator for +pbkdf
}

#[inline]
fn bitreverse(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

pub(crate) mod field {
    use super::bitreverse;
    use zkvm_jetpack::form::poly::Belt;

    #[inline]
    pub fn fft_in_place(data: &mut [Belt], log_n: usize) {
        let n = data.len() as u32;
        if n == 0 { return; }
        debug_assert_eq!(1u32 << log_n, n);

        for k in 0..n {
            let rk = bitreverse(k, log_n as u32);
            if k < rk { data.swap(k as usize, rk as usize); }
        }

        let order = Belt(n as u64);
        let root = order.ordered_root().unwrap_or(Belt::one());

        let mut m = 1u32;
        for _ in 0..log_n {
            let w_m = Belt::from(zkvm_jetpack::form::math::base::bpow(root.0, (n / (2 * m)) as u64));

            let mut k = 0;
            while k < n {
                let mut w = Belt::one();

                for j in 0..m {
                    let u = data[(k + j) as usize];
                    let v = data[(k + j + m) as usize] * w;
                    data[(k + j) as usize] = u + v;
                    data[(k + j + m) as usize] = u - v;
                    w = w * w_m;
                }

                k += 2 * m;
            }

            m *= 2;
        }
    }

    #[inline]
    pub fn ifft_in_place(data: &mut [Belt], log_n: usize) {
        let n = data.len() as u32;
        if n == 0 { return; }
        debug_assert_eq!(1u32 << log_n, n);

        for k in 0..n {
            let rk = bitreverse(k, log_n as u32);
            if k < rk { data.swap(k as usize, rk as usize); }
        }

        let order = Belt(n as u64);
        let inv_root = order.ordered_root().unwrap_or(Belt::one()).inv();

        let mut m = 1u32;
        for _ in 0..log_n {
            let w_m = Belt::from(zkvm_jetpack::form::math::base::bpow(inv_root.0, (n / (2 * m)) as u64));

            let mut k = 0;
            while k < n {
                let mut w = Belt::one();

                for j in 0..m {
                    let u = data[(k + j) as usize];
                    let v = data[(k + j + m) as usize] * w;
                    data[(k + j) as usize] = u + v;
                    data[(k + j + m) as usize] = u - v;
                    w = w * w_m;
                }

                k += 2 * m;
            }

            m *= 2;
        }

        let inv_n = order.inv();
        for x in data.iter_mut() {
            *x = *x * inv_n;
        }
    }

    pub use Belt;
}

#[cfg(test)]
mod tests {
    use super::{jet_stark_engine_jet_hook__compute_lde, field::{ifft_in_place, Belt}};

    #[test]
    fn test_fft_roundtrip() {
        let log_n = 3;
        let n = 1 << log_n;
        let mut buf: Vec<u64> = Vec::with_capacity(2 + n);
        buf.push(n as u64);
        buf.push(log_n as u64);
        for i in 0..n {
            buf.push((i as u64) + 1);
        }

        jet_stark_engine_jet_hook__compute_lde(buf.as_mut_ptr());

        let mut out: Vec<Belt> = buf[2..].iter().map(|&c| Belt(c)).collect();
        ifft_in_place(&mut out, log_n);
        let recovered: Vec<u64> = out.into_iter().map(|b| b.0).collect();

        let expected: Vec<u64> = (1..=n as u64).collect();
        assert_eq!(recovered, expected);
    }
}
