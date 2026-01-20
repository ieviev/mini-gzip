use std::{
    env,
    io::{self, Read, Write},
};

const MAXBITS: usize = 15;
const MAXLCODES: usize = 286;
const MAXDCODES: usize = 30;
const FIXLCODES: usize = 288;
const MAXDIST: usize = 32768;
const LENS: [u16; 29] = [
    3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131,
    163, 195, 227, 258,
];
const LEXT: [u16; 29] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0,
];
const DISTS: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537,
    2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
];
const DEXT: [u16; 30] = [
    0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13,
    13,
];
const ORDER: [usize; 19] = [
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15,
];

struct State<'a> {
    bit_count: i32,
    bit_buffer: i32,
    input: &'a [u8],
    pos: usize,
    next: usize,
    window: [u8; MAXDIST],
}

struct Huffman<'a> {
    count: &'a mut [i16],
    symbol: &'a mut [i16],
}

impl<'a> State<'a> {
    fn nextbyte(&mut self) -> u8 {
        if self.pos < self.input.len() {
            self.pos += 1;
            return self.input[self.pos - 1];
        }
        panic!("unexpected end of input");
    }

    fn bits(&mut self, need: i32) -> i32 {
        let mut val = self.bit_buffer;
        while self.bit_count < need {
            val |= (self.nextbyte() as i32) << self.bit_count;
            self.bit_count += 8;
        }
        self.bit_buffer = val >> need;
        self.bit_count -= need;
        val & ((1i32 << need) - 1)
    }

    fn output(&mut self, out: &mut Vec<u8>, c: u8) {
        out.push(c);
        self.window[self.next] = c;
        self.next = (self.next + 1) & (MAXDIST - 1);
    }
}

fn decode(s: &mut State, h: &Huffman) -> i32 {
    let (mut bitbuf, mut left) = (s.bit_buffer, s.bit_count);
    let (mut len, mut code, mut first, mut index) = (1, 0, 0, 0);
    let mut next_idx = 1usize;
    loop {
        while left > 0 {
            left -= 1;
            code |= bitbuf & 1;
            bitbuf >>= 1;
            let count = h.count[next_idx] as i32;
            next_idx += 1;
            if code < first + count {
                s.bit_buffer = bitbuf;
                s.bit_count = (s.bit_count - len) & 7;
                return h.symbol[(index + code - first) as usize] as i32;
            }
            index += count;
            first = (first + count) << 1;
            code <<= 1;
            len += 1;
        }
        left = (MAXBITS as i32 + 1) - len;
        if left == 0 {
            break;
        }
        bitbuf = s.nextbyte() as i32;
        if left > 8 {
            left = 8;
        }
    }
    panic!("invalid huffman code");
}

fn construct(h: &mut Huffman, length: &[u16], n: usize) -> i32 {
    h.count[..=MAXBITS].fill(0);
    for i in 0..n {
        h.count[length[i] as usize] += 1;
    }
    if h.count[0] as usize == n {
        return 0;
    }
    let mut left = 1i32;
    for len in 1..=MAXBITS {
        left = (left << 1) - h.count[len] as i32;
        if left < 0 {
            return left;
        }
    }
    let mut offs = [0i16; MAXBITS + 1];
    for len in 1..MAXBITS {
        offs[len + 1] = offs[len] + h.count[len];
    }
    for sym in 0..n {
        if length[sym] != 0 {
            h.symbol[offs[length[sym] as usize] as usize] = sym as i16;
            offs[length[sym] as usize] += 1;
        }
    }
    left
}

fn codes(s: &mut State, out: &mut Vec<u8>, lc: &Huffman, dc: &Huffman) {
    loop {
        let sym = decode(s, lc);
        if sym < 256 {
            s.output(out, sym as u8);
        } else if sym > 256 {
            let idx = (sym - 257) as usize;
            let len = LENS[idx] as i32 + s.bits(LEXT[idx] as i32);
            let dsym = decode(s, dc) as usize;
            let dist = DISTS[dsym] as u32 + s.bits(DEXT[dsym] as i32) as u32;
            for _ in 0..len {
                let c = s.window[s.next.wrapping_sub(dist as usize) & (MAXDIST - 1)];
                s.output(out, c);
            }
        } else {
            return;
        }
    }
}

fn fixed(s: &mut State, out: &mut Vec<u8>) {
    let (mut lc_cnt, mut lc_sym) = ([0i16; MAXBITS + 1], [0i16; FIXLCODES]);
    let (mut dc_cnt, mut dc_sym) = ([0i16; MAXBITS + 1], [0i16; MAXDCODES]);
    let mut lc = Huffman {
        count: &mut lc_cnt,
        symbol: &mut lc_sym,
    };
    let mut dc = Huffman {
        count: &mut dc_cnt,
        symbol: &mut dc_sym,
    };
    let mut lengths = [8u16; FIXLCODES];
    lengths[144..256].fill(9);
    lengths[256..280].fill(7);
    construct(&mut lc, &lengths, FIXLCODES);
    lengths[..MAXDCODES].fill(5);
    construct(&mut dc, &lengths, MAXDCODES);
    codes(s, out, &lc, &dc)
}

fn dynamic(s: &mut State, out: &mut Vec<u8>) {
    let (nlen, ndist, ncode) = (
        s.bits(5) as usize + 257,
        s.bits(5) as usize + 1,
        s.bits(4) as usize + 4,
    );
    let mut lengths = [0u16; MAXLCODES + MAXDCODES];
    for i in 0..ncode {
        lengths[ORDER[i]] = s.bits(3) as u16;
    }
    let (mut lc_cnt, mut lc_sym) = ([0i16; MAXBITS + 1], [0i16; MAXLCODES]);
    let mut lc = Huffman {
        count: &mut lc_cnt,
        symbol: &mut lc_sym,
    };
    construct(&mut lc, &lengths, 19);
    let mut idx = 0usize;
    while idx < nlen + ndist {
        let sym = decode(s, &lc);
        if sym < 16 {
            lengths[idx] = sym as u16;
            idx += 1;
        } else {
            let (len, rep) = match sym {
                16 => (lengths[idx - 1], 3 + s.bits(2) as usize),
                17 => (0, 3 + s.bits(3) as usize),
                _ => (0, 11 + s.bits(7) as usize),
            };
            for _ in 0..rep {
                lengths[idx] = len;
                idx += 1;
            }
        }
    }
    construct(&mut lc, &lengths, nlen);
    let (mut dc_cnt, mut dc_sym) = ([0i16; MAXBITS + 1], [0i16; MAXDCODES]);
    let mut dc = Huffman {
        count: &mut dc_cnt,
        symbol: &mut dc_sym,
    };
    construct(&mut dc, &lengths[nlen..], ndist);
    codes(s, out, &lc, &dc)
}

fn stored(s: &mut State, out: &mut Vec<u8>) {
    s.bits(s.bit_count);
    let len = s.bits(16) as u32;
    s.bits(16);
    for _ in 0..len {
        let c = s.nextbyte();
        s.output(out, c);
    }
}

pub fn inflate(input: &[u8]) -> Vec<u8> {
    let mut s = State {
        bit_count: 0,
        bit_buffer: 0,
        input,
        pos: 0,
        next: 0,
        window: [0; MAXDIST],
    };
    let mut out = Vec::new();
    loop {
        let last = s.bits(1);
        match s.bits(2) {
            0 => stored(&mut s, &mut out),
            1 => fixed(&mut s, &mut out),
            2 => dynamic(&mut s, &mut out),
            _ => panic!("invalid block type"),
        }
        if last != 0 {
            break;
        }
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let buf: Vec<u8> = if args.len() > 1 {
        std::fs::read(&args[1]).expect("read file")
    } else {
        let mut b = Vec::new();
        io::stdin().read_to_end(&mut b).expect("read stdin");
        b
    };
    assert!(buf[0] == 0x1F && buf[1] == 0x8B, "not a gzip file");
    let mut offset = 10;
    if buf[3] == 8 {
        offset += buf[10..].iter().position(|&b| b == 0).unwrap() + 1;
    } else if buf[3] != 0 {
        panic!("unsupported gzip format");
    }
    io::stdout().write_all(&inflate(&buf[offset..])).unwrap();
}
