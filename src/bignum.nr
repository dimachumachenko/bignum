use crate::limbs::{LIMB_BITS, normalize_limb, normalize_limbs};

pub struct BigNum<let N: u32> {
    pub limbs: [u128; N],
}

impl<let N: u32> BigNum<N> {
    pub fn zero() -> Self {
        BigNum { limbs: [0u128; N] }
    }

    pub fn from_limbs(l: [u128; N]) -> Self {
        BigNum { limbs: normalize_limbs::<N>(l) }
    }

    pub fn from_bytes_le(bytes: [u8]) -> Self {
        // Pack LE bytes into 120-bit limbs
        let mut out: [u128; N] = [0u128; N];
        let mut bitpos: u32 = 0;
        let mut limb_idx: u32 = 0;

        for i in 0..bytes.len() {
            let b128: u128 = (bytes[i] as u128);
            let shift: u32 = (bitpos % LIMB_BITS);

            if limb_idx < N {
    		out[limb_idx] = normalize_limb(out[limb_idx] | (b128 << (shift as u128)));
	    }


            bitpos += 8u32;

            if (bitpos % LIMB_BITS) == 0u32 {
                limb_idx += 1u32;
                if limb_idx < N {
                    out[limb_idx] = 0u128;
                }
            } else {
                if shift > (LIMB_BITS - 8u32) {
                    // carry bits into the next limb
                    let hi_shift: u32 = (shift + 8u32) - LIMB_BITS;
                    limb_idx += 1u32;
                    if limb_idx < N {
                        let carry_bits: u128 = b128 >> ((8u32 - hi_shift) as u128);
                        out[limb_idx] = normalize_limb(out[limb_idx] | carry_bits);
                    }
                }
            }
        }
        BigNum { limbs: out }
    }

    pub fn to_bytes_le<let OUT: u32>(self) -> [u8; OUT] {
        let mut out: [u8; OUT] = [0u8; OUT];

        let mut bitpos: u32 = 0;
        for i in 0..OUT {
            let limb_idx: u32 = bitpos / LIMB_BITS;
            let shift: u32    = bitpos % LIMB_BITS;
            let mut byte: u128 = 0u128;

            if limb_idx < N {
                byte = (self.limbs[limb_idx] >> (shift as u128)) & 0xFFu128;

                if shift > (LIMB_BITS - 8u32) {
                    if (limb_idx + 1u32) < N {
                        let hi_shift: u32 = (shift + 8u32) - LIMB_BITS;
                        let mask_hi: u128 = (1u128 << (hi_shift as u128)) - 1u128;
                        let next = self.limbs[limb_idx + 1u32] & mask_hi;
                        byte |= next << ((8u32 - hi_shift) as u128);
                    }
                }
            }

            out[i] = (byte & 0xFFu128) as u8;
            bitpos += 8u32;
        }
        out
    }

    pub fn ct_select(a: Self, b: Self, flag: u128) -> Self {
        // mask is all-zeros or all-ones
        let mask: u128 = (!0u128) * (flag & 1u128);
        let mut out_arr: [u128; N] = [0u128; N];
        for i in 0..N {
            let ai = a.limbs[i];
            let bi = b.limbs[i];
            // out = ai ^ ((ai ^ bi) & mask)
            out_arr[i] = ai ^ ((ai ^ bi) & mask);
        }
        BigNum { limbs: out_arr }
    }
}

