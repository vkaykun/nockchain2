#[no_mangle]
pub extern "C" fn jet_compute_codeword_commitments(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-codeword-commitments
}

#[no_mangle]
pub extern "C" fn jet_stark_engine_jet_hook__compute_lde(ptr: *mut u64) {
    use core::slice;
    use zkvm_jetpack::form::math::base::*;
    use zkvm_jetpack::form::poly::Belt;

    // SAFETY: caller ensures `ptr` is valid for `len + 2` u64 words.
    unsafe {
        let len = ptr.read() as usize;
        let log_n = ptr.add(1).read() as u32;
        let coeff_ptr = ptr.add(2);
        let coeffs = slice::from_raw_parts_mut(coeff_ptr, len);

        let mut tmp: Vec<Belt> = Vec::with_capacity(len);
        for i in 0..len {
            tmp.push(Belt(coeffs[i]));
        }

        debug_assert_eq!(len as u32, 1u32 << log_n);
        let root = match Belt(len as u64).ordered_root() {
            Ok(r) => r,
            Err(_) => return,
        };

        radix2_fft_in_place(&mut tmp, &root);

        for i in 0..len {
            coeffs[i] = tmp[i].0;
        }
    }
}

fn radix2_fft_in_place(x: &mut [Belt], root: &Belt) {
    use zkvm_jetpack::form::math::base::bpow;
    let n = x.len() as u32;
    if n <= 1 {
        return;
    }
    let log_n = n.ilog2();

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            x.swap(k as usize, rk as usize);
        }
    }

    let mut m = 1u32;
    for _ in 0..log_n {
        let w_m: Belt = bpow(root.0, (n / (2 * m)) as u64).into();
        let mut k = 0u32;
        while k < n {
            let mut w = Belt(1);
            for j in 0..m {
                let u = x[(k + j) as usize];
                let v = x[(k + j + m) as usize] * w;
                x[(k + j) as usize] = u + v;
                x[(k + j + m) as usize] = u - v;
                w = w * w_m;
            }
            k += 2 * m;
        }
        m *= 2;
    }
}

fn bitreverse(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

#[no_mangle]
pub extern "C" fn jet_compute_composition_poly(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-composition-poly
}

#[no_mangle]
pub extern "C" fn jet_compute_deep(_ptr: *mut u64) {
    // TODO: implement native accelerator for +compute-deep
}

#[no_mangle]
pub extern "C" fn jet_precompute_ntts(_ptr: *mut u64) {
    // TODO: implement native accelerator for +precompute-ntts
}

#[no_mangle]
pub extern "C" fn jet_build_compute_queue(_ptr: *mut u64) {
    // TODO: implement native accelerator for +build-compute-queue
}

#[no_mangle]
pub extern "C" fn jet_build_jute_list(_ptr: *mut u64) {
    // TODO: implement native accelerator for +build-jute-list
}

#[no_mangle]
pub extern "C" fn jet_indirect_to_bits(_ptr: *mut u64) {
    // TODO: implement native accelerator for +indirect-to-bits
}

#[no_mangle]
pub extern "C" fn jet_noun_get_zero_mults(_ptr: *mut u64) {
    // TODO: implement native accelerator for +noun-get-zero-mults
}

#[no_mangle]
pub extern "C" fn jet_pbkdf(_ptr: *mut u64) {
    // TODO: implement native accelerator for +pbkdf
}
