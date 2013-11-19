import binascii
import random
import sha
import sys
import math

def stringToLong(s):
    return long(binascii.hexlify(s), 16)

def longToString(n):
    s = "%x" % n
    if len(s) % 2 == 1:
        s = '0' + s
    return binascii.unhexlify(s)

def is_prime(n, rounds=50):
    """
    Miller-Rabin prime probability test
    """
    r = n - 1
    s = 0
    # need to represent n - 1 as 2^s * t (t is even)
    while (not r & 1): #r % 2 == 0
        r = r >> 1 #r = r / 2
        s = s + 1
    for _ in xrange(rounds):
        # calculating a^t mod m
        a = random.randint(2, n - 1)
        y = pow(a, r, n)
        if (y != 1 and y != n - 1):
            for _ in xrange(s - 1):
                if (y != n - 1):
                    y = y * y % n
                    if (y == 1):
                        return False #complex
            if (y != n-1):
                return False #complex
    return True #probably prime


def rng(k):
    try:
        return random.getrandbits(k)
    except Exception, e:
        print 'Error: Cannot get random number of %d bits. \nDetails: %s' % (k, e.message)
        sys.exit(1)         


def sha1(x):
    s = longToString(x)
    return stringToLong(sha.new(s).digest())

def print_and_exit(p, q):
    print 'P:', p
    print 'Q:', q
    print 'P is prime:', is_prime(p)
    print 'Q is prime:', is_prime(q)
    sys.exit(0)

def weighted_sum(v, N, b):
    """
    Calculates v0 + v1 * 2^N + ... + v(n-1) * 2^(N*(n-1)) + 
    (vn mod 2b) * 2^(n*N)
    where N comes from desired numbers' restrictions
    and b comes from N's representation
    """
    w = v[0]
    n = len(v)- 1 # v consists of n+1 items, where n comes from L representation
    for i in xrange(1, n - 1):
        w += v[i] * (2 ** (N * i))
    divider = 2 ** (2 * b)
    w += (v[n] % divider) * (2 ** (n * N))
    return w

def calculate_seedlen(N):
    """
    Seed length should be >= N, and it will be better if it's a degree of 2
    So returning nearest degree of 2 that's more than N
    """
    deg = math.ceil(math.log(N, 2))
    intdeg = int(deg)
    return 2 ** intdeg

def helper_sequence(seed, offset, n):
    seedlen = seed.bit_length()

    v = [0] * (n+1)
    for k in range(0, n + 1):
        v[k] = sha1((seed + offset + k) % (2 ** seedlen))
    return v

def q_candidate(seed, N):
    """
    Calculate probable least of two desired numbers p and q - q (depends on N)
    (needs primality checking)
    """    
    seedlen = seed.bit_length()
    u = sha1(seed) ^ sha1((seed + 1) % (2 ** seedlen))
    # bitwise setting first and last bits to 1
    q = u | (2 ** (N - 1)) | 1
    return q

def seed_and_q(seedlen, N):
    """
    Calculates q and seed (bit length seedlen) with which q is probably prime
    q depends on N
    """
    while True:
        seed = rng(seedlen)
        q = q_candidate(seed, N)   
        if is_prime(q):
            return (seed, q)

def p_candidate(L, N, seed, offset, n, b):
    v = helper_sequence(seed, offset, n)
        
    w = weighted_sum(v, N, b)
    x = w + 2 ** (L - 1)
    c = x % (2 * q)
    p = x - (c - 1) 
    return p

# rules:
# 2^(n-1) < q < 2^n
# 2^(L-1) < p < 2^l
# p - 1 mod q == 0
if len(sys.argv) < 3:
    print 'Usage: ./fips_primes.py <L> <N>'
    print 'Generates two possibly prime numbers p and q with given restrictions : '
    print ' <N>                 - 2 ^ (N - 1) < q < 2 ^ N'
    print ' <L>                 - 2 ^ (L - 1) < p < 2 ^ L'
    print ' p - 1 mod q == 0 '
    print ' Example input pairs: L = 1024 N = 160 '
    print '                      L = 2048 N = 224 '
    print '                      L = 2048 N = 256 '
    print ' N should be > 4. L should be > 20'
    sys.exit(1)

try:
    N = int(sys.argv[2])
    if N <= 4:
        raise
    print 'Got N = %d' % N
except:
    print 'Invalid N'
    sys.exit(1)

try:
    L = int(sys.argv[1])
    if L <= 20:
        raise
    print 'Got L = %d' % L
except:
    print 'Invalid L'
    sys.exit(1)

if N > L:
    print 'N should be less than L'
    sys.exit(1)
# L = 1024
# N = 64

# let L be: L - 1 == n * N + b, where 0 <= b < N
(n, b) = divmod(L, N)
seedlen = calculate_seedlen(N)

# while True:
for counter in xrange(0, 2 ** 4): # stop after 2 ** 12 loops if can't get primes
    seed, q = seed_and_q(seedlen, N)     

    offset = 2 # magic number from algo

    while True:
        p = p_candidate(L, N, seed, offset, n, b)

        if p >= (2 ** (L - 1)):
            if is_prime(p):
                print_and_exit(p, q)

        offset += n + 1


