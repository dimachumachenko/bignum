use crate::limbs::{LIMB_BITS, MASK_120, normalize_limb, adc_120, sbb_120, add_no_mod, sub_no_mod, cmp_limbs, normalize_limbs};

#[test]
fn globals_and_mask_sanity() {
    println("== Global sanity check ==");
    println("LIMB_BITS =");
    println(LIMB_BITS);

    println("MASK_120 =");
    println(MASK_120);

    assert(LIMB_BITS == 120);
    assert(((MASK_120 + 1u128) & MASK_120) == 0u128);
    assert(normalize_limb(MASK_120) == MASK_120);

    let with_garbage = MASK_120 | (1u128 << 127u128);
    println("normalize_limb(garbage) =");
    println(normalize_limb(with_garbage));

    assert(normalize_limb(with_garbage) == MASK_120);
}

#[test]
fn adc_edges() {
    println("== adc_edges ==");

    let (lo1, c1) = adc_120(5u128, 7u128, 0u128);
    println("case1: lo");
    println(lo1);
    println("case1: c");
    println(c1);
    assert(lo1 == 12u128);
    assert(c1 == 0u128);

    let (lo2, c2) = adc_120(MASK_120, 1u128, 0u128);
    println("case2: lo");
    println(lo2);
    println("case2: c");
    println(c2);
    assert(lo2 == 0u128);
    assert(c2 == 1u128);

    let (lo3, c3) = adc_120(MASK_120, MASK_120, 1u128);
    println("case3: lo");
    println(lo3);
    println("case3: c");
    println(c3);
    assert(lo3 == MASK_120);
    assert(c3 == 1u128);

    let garbage_a = MASK_120 | (1u128 << 127u128);
    let garbage_b = MASK_120 | (1u128 << 126u128);
    let (lo4, c4) = adc_120(garbage_a, garbage_b, 0u128);
    println("case4: lo");
    println(lo4);
    println("case4: c");
    println(c4);
    assert(lo4 == MASK_120 - 1u128);
    assert(c4 == 1u128);
}

#[test]
fn adc_invariant_reconstruction() {
    println("== adc_invariant_reconstruction ==");

    let a = MASK_120;
    let b = (1u128 << 60u128);
    let cin = 1u128;

    let (lo, c) = adc_120(a, b, cin);

    println("lo");
    println(lo);
    println("carry");
    println(c);

    let left  = lo + (c << (LIMB_BITS as u128));
    let right = (a & MASK_120) + (b & MASK_120) + (cin & 1u128);

    println("left");
    println(left);
    println("right");
    println(right);

    assert(left == right);
}


#[test]
fn sbb_edges() {
    println("== sbb_edges ==");

    let (d1, b1) = sbb_120(5u128, 3u128, 0u128);
    println("case1 d"); 
    println(d1); 
    println("case1 borrow"); 
    println(b1);
    assert(d1 == 2u128);
    assert(b1 == 0u128);

    let (d2, b2) = sbb_120(0u128, 1u128, 0u128);
    println("case2 d"); 
    println(d2); 
    println("case2 borrow"); 
    println(b2);
    assert(d2 == ((1u128 << (LIMB_BITS as u128)) - 1u128));
    assert(b2 == 1u128);

    let (d3, b3) = sbb_120(10u128, 10u128, 1u128);
    println("case3 d"); 
    println(d3); 
    println("case3 borrow"); 
    println(b3);
    assert(d3 == ((1u128 << (LIMB_BITS as u128)) - 1u128));
    assert(b3 == 1u128);
}

#[test]
fn add_sub_multi_limb_roundtrip() {
    println("== add_sub_multi_limb_roundtrip ==");

    let a: [u128; 2] = [MASK_120, 0];
    let b: [u128; 2] = [1u128, 0];
    let (sum, c) = add_no_mod::<2>(a, b);
    println("sum[0]"); 
    println(sum[0]); 
    println("sum[1]"); 
    println(sum[1]); 
    println("carry"); 
    println(c);
    assert(sum[0] == 0u128);
    assert(sum[1] == 1u128);
    assert(c == 0u128);

    let (back, bor) = sub_no_mod::<2>(sum, b);
    println("back[0]"); 
    println(back[0]); 
    println("back[1]"); 
    println(back[1]); 
    println("borrow"); 
    println(bor);
    assert(bor == 0u128);
    assert(back[0] == a[0]);
    assert(back[1] == a[1]);
}

#[test]
fn sub_chain_borrow_across_limbs() {
    println("== sub_chain_borrow_across_limbs ==");

    let a: [u128; 2] = [0u128, 1u128];
    let b: [u128; 2] = [1u128, 0u128];
    let (d, br) = sub_no_mod::<2>(a, b);
    println("d[0]"); 
    println(d[0]); 
    println("d[1]"); 
    println(d[1]); 
    println("borrow"); 
    println(br);

    let LIMB_BASE = (1u128 << (LIMB_BITS as u128));
    assert(d[0] == (LIMB_BASE - 1u128));
    assert(d[1] == 0u128);
    assert(br == 0u128);
}

#[test]
fn cmp_and_normalize() {
    println("== cmp_and_normalize ==");

    let a: [u128; 3] = [0, 0, 1];
    let b: [u128; 3] = [0, 0, 1];
    let c: [u128; 3] = [0, 0, 2];
    let d: [u128; 3] = [1, 0, 0];

    println("cmp a,b"); 
    println(cmp_limbs::<3>(a, b));
    assert(cmp_limbs::<3>(a, b) == 0);
    assert(cmp_limbs::<3>(a, c) == -1);
    assert(cmp_limbs::<3>(c, a) == 1);
    assert(cmp_limbs::<3>(d, a) == -1); 

    let x: [u128; 2] = [MASK_120 | (1u128 << 127u128), MASK_120 | (1u128 << 126u128)];
    let y = normalize_limbs::<2>(x);
    assert(y[0] == MASK_120);
    assert(y[1] == MASK_120);
}

