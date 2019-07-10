uint uhash1(uint x) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}

uint uhash2(uvec2 v) { return uhash1(v.x ^ uhash1(v.y)); }
uint uhash3(uvec3 v) { return uhash1(v.x ^ uhash1(v.y ^ uhash1(v.z))); }
uint uhash4(uvec4 v) { return uhash1(v.x ^ uhash1(v.y ^ uhash1(v.z ^ uhash1(v.w)))); }

uint ihash1(int x)   { return uhash1(uint(x)); }
uint ihash2(ivec2 v) { return uhash2(uvec2(v)); }
uint ihash3(ivec3 v) { return uhash3(uvec3(v)); }
uint ihash4(ivec4 v) { return uhash4(uvec4(v)); }

uint hash1(float x) { return uhash1(floatBitsToUint(x)); }
uint hash2(vec2 v)  { return uhash2(floatBitsToUint(v)); }
uint hash3(vec3 v)  { return uhash3(floatBitsToUint(v)); }
uint hash4(vec4 v)  { return uhash4(floatBitsToUint(v)); }

// Construct a float with half-open range [0:1] using low 23 bits.
// All zeroes yields 0.0, all ones yields the next smallest representable value below 1.0.
float floatConstructUnsigned(uint m) {
    const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    m |= ieeeOne;                          // Add fractional part to 1.0

    float  f = uintBitsToFloat( m );       // Range [1:2]
    return f - 1.0;                        // Range [0:1]
}

float floatConstructSigned(uint m) {
    return floatConstructUnsigned(m) * 2.0 - 1.0;
}

// These function return a deterministically random floating point value in the
// range (0; 1)
float uhash1AsUnsignedFloat(uint x)  { return floatConstructUnsigned(uhash1(x)); }
float uhash2AsUnsignedFloat(uvec2 v) { return floatConstructUnsigned(uhash2(v)); }
float uhash3AsUnsignedFloat(uvec3 v) { return floatConstructUnsigned(uhash3(v)); }
float uhash4AsUnsignedFloat(uvec4 v) { return floatConstructUnsigned(uhash4(v)); }

float ihash1AsUnsignedFloat(int x)   { return floatConstructUnsigned(ihash1(x)); }
float ihash2AsUnsignedFloat(ivec2 v) { return floatConstructUnsigned(ihash2(v)); }
float ihash3AsUnsignedFloat(ivec3 v) { return floatConstructUnsigned(ihash3(v)); }
float ihash4AsUnsignedFloat(ivec4 v) { return floatConstructUnsigned(ihash4(v)); }

float hash1AsUnsignedFloat(float x) { return floatConstructUnsigned(hash1(x)); }
float hash2AsUnsignedFloat(vec2 v)  { return floatConstructUnsigned(hash2(v)); }
float hash3AsUnsignedFloat(vec3 v)  { return floatConstructUnsigned(hash3(v)); }
float hash4AsUnsignedFloat(vec4 v)  { return floatConstructUnsigned(hash4(v)); }

// These function return a deterministically random floating point value in the
// range (-1; 1)
float uhash1AsSignedFloat(uint x)  { return floatConstructSigned(uhash1(x)); }
float uhash2AsSignedFloat(uvec2 v) { return floatConstructSigned(uhash2(v)); }
float uhash3AsSignedFloat(uvec3 v) { return floatConstructSigned(uhash3(v)); }
float uhash4AsSignedFloat(uvec4 v) { return floatConstructSigned(uhash4(v)); }

float ihash1AsSignedFloat(int x)   { return floatConstructSigned(ihash1(x)); }
float ihash2AsSignedFloat(ivec2 v) { return floatConstructSigned(ihash2(v)); }
float ihash3AsSignedFloat(ivec3 v) { return floatConstructSigned(ihash3(v)); }
float ihash4AsSignedFloat(ivec4 v) { return floatConstructSigned(ihash4(v)); }

float hash1AsSignedFloat(float x) { return floatConstructSigned(hash1(x)); }
float hash2AsSignedFloat(vec2 v)  { return floatConstructSigned(hash2(v)); }
float hash3AsSignedFloat(vec3 v)  { return floatConstructSigned(hash3(v)); }
float hash4AsSignedFloat(vec4 v)  { return floatConstructSigned(hash4(v)); }

/*
 * Implementations of Quantile (inverse CDF) functions
 */

vec2 boxMuller2to2(vec2 uniformVector) {
    float term1 = sqrt(-2 * log(uniformVector.x));
    float term2 = TAU * uniformVector.y;

    return vec2(
        term1 * cos(term2),
        term1 * sin(term2)
    );
}

float boxMuller2to1(vec2 uniformVector) {
    float term1 = sqrt(-2 * log(uniformVector.x));
    float term2 = TAU * uniformVector.y;

    return term1 * cos(term2);
}

// An implementation of pseudo random normal variable sampling using the
// box-muller method
float uboxMullerSample(uint seed) {
    uint xHash = uhash1(seed);
    uint yHash = uhash1(xHash);
    float x = floatConstructUnsigned(xHash);
    float y = floatConstructUnsigned(yHash);

    return boxMuller2to1(vec2(x, y));
}

float iboxMullerSample(int seed) { return uboxMullerSample(uint(seed)); }
float boxMullerSample(float seed) { return uboxMullerSample(floatBitsToUint(seed)); }

/*
 * These two following methods `branchlessNormalQuantileApprox` and
 * `hybridNormalQuantileApprox` were taken from _Fast and Accurate Parallel
 * Computation of Quantile Functions for Random Number Generation_ by _Thomas
 * Luu_.
 */
float branchlessNormalQuantileApprox(float u) {
    float ushift = u - 0.5f;
    if (ushift > 0.0f) u = 1.0f - u;
    float v = -log(u + u);
    float p =   1.68267776058639e-6f;
    p = p * v + 0.0007404314351202936f;
    p = p * v + 0.03602364419560667f;
    p = p * v + 0.4500443083534446f;
    p = p * v + 1.861100468283588f;
    p = p * v + 2.748475794390544f;
    p = p * v + 1.253314132218524f;
    float q =   0.00003709787159774307f;
    q = q * v + 0.004513659269519104f;
    q = q * v + 0.1101701640048184f;
    q = q * v + 0.8410203004476538f;
    q = q * v + 2.402969434512837f;
    q = q * v + 2.692965915568952f;
    q = q * v + 1.0f;
    /* return __fdividef(p, q) * copysignf(v, ushift); */
    return (p / q) * (abs(v) * sign(ushift));
}

float hybridNormalQuantileApprox(float u) {
    float v, p, q, ushift;
    ushift = u - 0.5f;
    /* v = copysignf(ushift, 0.0f); */
    v = abs(ushift);
    if (v < 0.499433f) {
        /* asm("rsqrt.approx.ftz.f32 %0,%1;" : "=f"(v) : "f"(u - u * u)); */
        v = inversesqrt(u - u * u);
        v *= 0.5f;
        p =         0.001732781974270904f;
        p = p * v + 0.1788417306083325f;
        p = p * v + 2.804338363421083f;
        p = p * v + 9.35716893191325f;
        p = p * v + 5.283080058166861f;
        p = p * v + 0.07885390444279965f;
        p *= ushift;
        q =         0.0001796248328874524f;
        q = q * v + 0.02398533988976253f;
        q = q * v + 0.4893072798067982f;
        q = q * v + 2.406460595830034f;
        q = q * v + 3.142947488363618f;
    } else {
        if (ushift > 0.0f) u = 1.0f - u;
        /* asm("lg2.approx.ftz.f32 %0,%1;" : "=f"(v)) : "f"(u + u)); */
        v = log2(u + u);
        v *= -0.6931471805599453f;
        if (v < 22.0f) {
            p =         0.000382438382914666f;
            p = p * v + 0.03679041341785685f;
            p = p * v + 0.5242351532484291f;
            p = p * v + 1.21642047402659f;
            q =         9.14019972725528e-6f;
            q = q * v + 0.003523083799369908f;
            q = q * v + 0.126802543865968f;
            q = q * v + 0.8502031783957995f;
        } else {
            p =         0.00001016962895771568f;
            p = p * v + 0.003330096951634844f;
            p = p * v + 0.1540146885433827f;
            p = p * v + 1.045480394868638f;
            q =         1.303450553973082e-7f;
            q = q * v + 0.0001728926914526662f;
            q = q * v + 0.02031866871146244f;
            q = q * v + 0.3977137974626933f;
        }
        /* p *= copysignf(v, ushift); */
        p *= abs(v) * sign(ushift);
    }
    q = q * v + 1.0f;
    /* asm("rcp.approx.ftz.f32 %0,%1;" : "=f"(v) : "f"(q)); */
    v = 1 / q;
    return p * v;
}

float ubranchlessNormalQuantileApproxSample(uint seed) {
    float uniformVariable = uhash1AsUnsignedFloat(seed);
    return branchlessNormalQuantileApprox(uniformVariable);
}
float ibranchlessNormalQuantileApproxSample(int seed) { return ubranchlessNormalQuantileApproxSample(uint(seed)); }
float branchlessNormalQuantileApproxSample(float seed) { return ubranchlessNormalQuantileApproxSample(floatBitsToUint(seed)); }

float uhybridNormalQuantileApproxSample(uint seed) {
    float uniformVariable = uhash1AsUnsignedFloat(seed);
    return hybridNormalQuantileApprox(uniformVariable);
}
float ihybridNormalQuantileApproxSample(int seed) { return uhybridNormalQuantileApproxSample(uint(seed)); }
float hybridNormalQuantileApproxSample(float seed) { return uhybridNormalQuantileApproxSample(floatBitsToUint(seed)); }

/*
 * Shorthands for the preferred way of sampling normally distributed random
 * variables
 */
float unormal(uint seed) { return uboxMullerSample(seed); }
float inormal(int seed) { return iboxMullerSample(seed); }
float normal(float seed) { return boxMullerSample(seed); }

/*
 * Methods to generate a random unit vector. Implementations up to 3 dimensions
 * use the closed-form formula, which does not exist for higher dimensions,
 * where we have to use normally distributed random variable sampling.
 */
float ihash1AsUnitVec1(int x) {
    uint hash = ihash1(x);
    return ((hash & 1) == 0) ? -1.0 : 1.0;
}

vec2 ihash2AsUnitVec2(ivec2 v) {
    float hash = ihash2AsUnsignedFloat(v);
    float angle = hash * TAU;
    return vec2(cos(angle), sin(angle));
}

vec3 ihash3AsUnitVec3(ivec3 v) {
    uint hash1 = ihash3(v);
    uint hash2 = uhash1(hash1);
    float phi = floatConstructUnsigned(hash1) * TAU;
    float theta = acos(2 * floatConstructUnsigned(hash2) - 1);
    return vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}

vec4 ihash4AsUnitVec4(ivec4 v) {
    const float maxCoord = uintBitsToFloat(0x7f7fffff);
    const float minCoord = uintBitsToFloat(0xff7fffff);
    uint hash1 = ihash4(v);
    uint hash2 = uhash1(hash1);
    uint hash3 = uhash1(hash2);
    uint hash4 = uhash1(hash3);
    float normal1 = unormal(hash1);
    float normal2 = unormal(hash2);
    float normal3 = unormal(hash3);
    float normal4 = unormal(hash4);
    vec4 vector = vec4(
        clamp(normal1, minCoord, maxCoord),
        clamp(normal2, minCoord, maxCoord),
        clamp(normal3, minCoord, maxCoord),
        clamp(normal4, minCoord, maxCoord)
    );

    return normalize(vector);
}

/*
 * Derivation:
 * mix(a, b, t) = a * (1 - t) + b * t
 * mix(a, b, t) = a - at + b * t
 *
 * dmix(a, b, t)/dt = d(a - at + b * t)/dt
 * dmix(a, b, t)/dt = -a + b
 *
 * dmix(a, b, t)/da = d(a - at + b * t)/da
 * dmix(a, b, t)/da = 1 - t
 *
 * dmix(a, b, t)/db = d(a - at + b * t)/db
 * dmix(a, b, t)/db = t
 */
float mixDerivativeT(float a, float b, float t) { return b - a; }
float mixDerivativeA(float a, float b, float t) { return 1 - t; }
float mixDerivativeB(float a, float b, float t) { return     t; }

// Value noise assigns random values to lattice points and interpolates between them
float valueNoiseUnsigned1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnsignedFloat(l0);
    float l1Hash = ihash1AsUnsignedFloat(l1);
    // Interpolation
    float t = fract(x);
    float tMapped = noiseInterpolation(t);

    return mix(l0Hash, l1Hash, tMapped);
}

/*
 * Hand-expanded derivation, see `valueNoiseUnsignedDerivative1` for a general
 * solution. See `smoothStepDerivative1` and `mixDerivativeT` for
 * solutions to partial derivatives of used functions.
 *
 * f(t) = mix(l0, l1, smoothStep(t))
 * f'(t) = mix'(l0, l1, smoothStep(t)) * smoothStep'(t)
 * f'(t) = (l1 - l0) * 6t (1 - t);
 */
float valueNoiseUnsignedDerivativeCustom1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnsignedFloat(l0);
    float l1Hash = ihash1AsUnsignedFloat(l1);
    // Interpolation
    float t = fract(x);

    return (l1Hash - l0Hash) * 6 * t * (1 - t);
}

float valueNoiseUnsignedGradient1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnsignedFloat(l0);
    float l1Hash = ihash1AsUnsignedFloat(l1);
    // Interpolation
    float t = fract(x);
    float tMapped = noiseInterpolation(t);
    float tMappedDerivative = noiseInterpolationGradient(t);

    /*
     * z = dmix(a, b, t2)
     * dmix(l0, l1, smoothStep(t))/dt =
     *   (dmix(l0, l1, smoothStep(t))/da) * (da/dt) (= 0)
     *   (dmix(l0, l1, smoothStep(t))/db) * (db/dt) (= 0)
     *   (dmix(l0, l1, smoothStep(t))/dt2) * (dt2/dt)
     * because
     *   da/dt = 0
     *   db/dt = 0
     */
    return mixDerivativeT(l0Hash, l1Hash, tMapped) * tMappedDerivative;
}

float valueNoiseUnsigned2(vec2 v) {
    // Lattice point integer coordinates
    ivec2 l00 = ivec2(floor(v));
    ivec2 l01 = l00 + ivec2(0, 1);
    ivec2 l10 = l00 + ivec2(1, 0);
    ivec2 l11 = l00 + ivec2(1, 1);
    // Deterministically random values for lattice points
    float l00Hash = ihash2AsUnsignedFloat(l00);
    float l01Hash = ihash2AsUnsignedFloat(l01);
    float l10Hash = ihash2AsUnsignedFloat(l10);
    float l11Hash = ihash2AsUnsignedFloat(l11);
    // Interpolation
    vec2 t = fract(v);
    vec2 tMapped = noiseInterpolation(t);

    return mix(
        mix(l00Hash, l01Hash, tMapped.y),
        mix(l10Hash, l11Hash, tMapped.y),
        tMapped.x
    );
}

vec2 valueNoiseUnsignedGradient2(vec2 v) {
    // Lattice point integer coordinates
    ivec2 l00 = ivec2(floor(v));
    ivec2 l01 = l00 + ivec2(0, 1);
    ivec2 l10 = l00 + ivec2(1, 0);
    ivec2 l11 = l00 + ivec2(1, 1);
    // Deterministically random values for lattice points
    float l00Hash = ihash2AsUnsignedFloat(l00);
    float l01Hash = ihash2AsUnsignedFloat(l01);
    float l10Hash = ihash2AsUnsignedFloat(l10);
    float l11Hash = ihash2AsUnsignedFloat(l11);
    // Interpolation
    vec2 t = fract(v);
    vec2 tMapped = noiseInterpolation(t);
    vec2 tMappedDerivative = noiseInterpolationGradient(t);

    /*
     * Differentiate with respect to t.x and t.y:
     *
     * z = dmix(a, b, t)
     *
     * dz/dtx:
     * dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/dtx =
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/da) * (da/dtx)
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/db) * (db/dtx)
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/dt) * (dt/dtx)
     * where
     * da/dtx = dmix(l00, l01, smoothStep(t).y)/dtx =
     * + (dmix(l00, l01, smoothStep(t).y)/dai) * (dai/dtx) (= 0)
     * + (dmix(l00, l01, smoothStep(t).y)/dbi) * (dbi/dtx) (= 0)
     * + (dmix(l00, l01, smoothStep(t).y)/dti) * (dti/dtx) (= 0)
     * db/dtx = analogously (= 0)
     * dt/dtx = d(smoothStep(t).x)/dtx
     *
     * dz/dty:
     * dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/dtx =
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/da) * (da/dty)
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/db) * (db/dty)
     * + (dmix(mix(l00, l01, smoothStep(t).y), mix(l00, l01, smoothStep(t).y), smoothStep(t).x)/dt) * (dt/dty)
     * where
     * da/dty = dmix(l00, l01, smoothStep(t).y)/dtx =
     * + (dmix(l00, l01, smoothStep(t).y)/dai) * (dai/dty) (= 0)
     * + (dmix(l00, l01, smoothStep(t).y)/dbi) * (dbi/dty) (= 0)
     * + (dmix(l00, l01, smoothStep(t).y)/dti) * (dti/dty)
     *      where
     *      dti/dty = d(smoothStep(t).y)/dty
     *
     * db/dty = analogously
     * dt/dty = d(smoothStep(t).x)/dty (= 0)
     */

    /* Direct rewrite without the usage of temporary variables:
    return vec2(
        mixDerivativeT(
            mix(l00Hash, l01Hash, tMapped.y),
            mix(l10Hash, l11Hash, tMapped.y),
            tMapped.x
        ) * tMappedDerivative.x,

        mixDerivativeA(
            mix(l00Hash, l01Hash, tMapped.y),
            mix(l10Hash, l11Hash, tMapped.y),
            tMapped.x
        ) * mixDerivativeT(l00Hash, l01Hash, tMapped.y) * tMappedDerivative.y
        + mixDerivativeB(
            mix(l00Hash, l01Hash, tMapped.y),
            mix(l10Hash, l11Hash, tMapped.y),
            tMapped.x
        ) * mixDerivativeT(l10Hash, l11Hash, tMapped.y) * tMappedDerivative.y
    );
    */

    float mixl0_Hash = mix(l00Hash, l01Hash, tMapped.y);
    float mixl1_Hash = mix(l10Hash, l11Hash, tMapped.y);

    return vec2(
        // dz/dtx
        mixDerivativeT(
            mixl0_Hash,
            mixl1_Hash,
            tMapped.x
        ) * tMappedDerivative.x,

        // dz/dty
        (
            mixDerivativeA(
                mixl0_Hash,
                mixl1_Hash,
                tMapped.x
            ) * mixDerivativeT(
                l00Hash,
                l01Hash,
                tMapped.y
            ) + mixDerivativeB(
                mixl0_Hash,
                mixl1_Hash,
                tMapped.x
            ) * mixDerivativeT(
                l10Hash,
                l11Hash,
                tMapped.y
            )
        ) * tMappedDerivative.y
    );
}

float valueNoiseUnsigned3(vec3 v) {
    // Lattice point integer coordinates
    ivec3 l000 = ivec3(floor(v));
    ivec3 l001 = l000 + ivec3(0, 0, 1);
    ivec3 l010 = l000 + ivec3(0, 1, 0);
    ivec3 l011 = l000 + ivec3(0, 1, 1);
    ivec3 l100 = l000 + ivec3(1, 0, 0);
    ivec3 l101 = l000 + ivec3(1, 0, 1);
    ivec3 l110 = l000 + ivec3(1, 1, 0);
    ivec3 l111 = l000 + ivec3(1, 1, 1);
    // Deterministically random values for lattice points
    float l000Hash = ihash3AsUnsignedFloat(l000);
    float l001Hash = ihash3AsUnsignedFloat(l001);
    float l010Hash = ihash3AsUnsignedFloat(l010);
    float l011Hash = ihash3AsUnsignedFloat(l011);
    float l100Hash = ihash3AsUnsignedFloat(l100);
    float l101Hash = ihash3AsUnsignedFloat(l101);
    float l110Hash = ihash3AsUnsignedFloat(l110);
    float l111Hash = ihash3AsUnsignedFloat(l111);
    // Interpolation
    vec3 t = fract(v);
    vec3 tMapped = noiseInterpolation(t);

    return mix(
        mix(
            mix(l000Hash, l001Hash, tMapped.z),
            mix(l010Hash, l011Hash, tMapped.z),
            tMapped.y
        ),
        mix(
            mix(l100Hash, l101Hash, tMapped.z),
            mix(l110Hash, l111Hash, tMapped.z),
            tMapped.y
        ),
        tMapped.x
    );
}

vec3 valueNoiseUnsignedGradient3(vec3 v) {
    // Lattice point integer coordinates
    ivec3 l000 = ivec3(floor(v));
    ivec3 l001 = l000 + ivec3(0, 0, 1);
    ivec3 l010 = l000 + ivec3(0, 1, 0);
    ivec3 l011 = l000 + ivec3(0, 1, 1);
    ivec3 l100 = l000 + ivec3(1, 0, 0);
    ivec3 l101 = l000 + ivec3(1, 0, 1);
    ivec3 l110 = l000 + ivec3(1, 1, 0);
    ivec3 l111 = l000 + ivec3(1, 1, 1);
    // Deterministically random values for lattice points
    float l000Hash = ihash3AsUnsignedFloat(l000);
    float l001Hash = ihash3AsUnsignedFloat(l001);
    float l010Hash = ihash3AsUnsignedFloat(l010);
    float l011Hash = ihash3AsUnsignedFloat(l011);
    float l100Hash = ihash3AsUnsignedFloat(l100);
    float l101Hash = ihash3AsUnsignedFloat(l101);
    float l110Hash = ihash3AsUnsignedFloat(l110);
    float l111Hash = ihash3AsUnsignedFloat(l111);
    // Interpolation
    vec3 t = fract(v);
    vec3 tMapped = noiseInterpolation(t);
    vec3 tMappedDerivative = noiseInterpolationGradient(t);

    /*
     * Differentiate with respect to t.x, t.y and t.z:
     *
     * z = dmix(a, b, t)
     *
     * dz/dtx:
     * dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dtx =
     *   (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/da) * (da/dtx) (= 0, see below)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/db) * (db/dtx) (= 0, see below)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dt) * (dt/dtx)
     *   where
     *   da/dtx = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dtx) =
     *     (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dai/dtx) (= 0, see below)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dbi/dtx) (= 0, see below)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (dti/dtx) (= 0, see below)
     *     where
     *     dai/dtx = dmix(l000, l001, smoothStep(t).z)/dtx = 0
     *     dbi/dtx = 0 (analogously)
     *     dti/dtx = 0 (analogously)
     *     therefore
     *     da/dtx = 0
     *   db/dtx = 0 (analogously)
     *   dt/dtx = d(smoothStep(t).x)/dtx
     *   therefore
     *   dz/dtx = (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dt) * (d(smoothStep(t).x)/dtx)
     *
     * dz/dty:
     * dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dty =
     *   (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/da) * (da/dty)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/db) * (db/dty)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dt) * (dt/dty)
     *   where
     *   da/dty = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dty) =
     *     (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dai/dty) (= 0, see below)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dbi/dty) (= 0, see below)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (dti/dty)
     *     where
     *     dai/dty = dmix(l000, l001, smoothStep(t).z)/dty = 0
     *     dbi/dty = 0 (analogously)
     *     dti/dty = d(smoothStep(t).y)/dty
     *     therefore
     *     da/dty = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (dti/dty)
     *            = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (d(smoothStep(t).y)/dty)
     *   db/dty = (analogously) = (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dti) * (dti/dty)
     *          = (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dti) * (d(smoothStep(t).y)/dty)
     *   dt/dty = d(smoothStep(t).x)/dty = 0
     *   therefore
     *   dz/dty =
     *     (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/da)
     *   * (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (d(smoothStep(t).y)/dty)
     *   + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/db)
     *   * (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dti) * (d(smoothStep(t).y)/dty)
     *
     * dz/dtz:
     * dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dtz =
     *   (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/da) * (da/dtz)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/db) * (db/dtz)
     * + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/dt) * (dt/dtz)
     *   where
     *   da/dtz = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dtz) =
     *     (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dai/dtz)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dbi/dtz)
     *   + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dti) * (dti/dtz) (= 0, see below)
     *     where
     *     dai/dtz = dmix(l000, l001, smoothStep(t).z)/dtz = (dmix(l000, l001, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)
     *     dbi/dtz = dmix(l010, l011, smoothStep(t).z)/dtz = (dmix(l010, l011, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)
     *     dti/dtz = d(smoothStep(t).y)/dtz = 0
     *     therefore
     *     da/dtz = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dai/dtz)
     *            + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dbi/dtz)
     *            = (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai)
     *            * (dmix(l000, l001, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)
     *            + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi)
     *            * (dmix(l010, l011, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz))
     *   db/dtz = (analogously) =
     *          = (dmix(mix(l100, l001, smoothStep(t).z), mix(l110, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dai/dtz)
     *          + (dmix(mix(l100, l001, smoothStep(t).z), mix(l110, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dbi/dtz)
     *          = (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dai)
     *          * (dmix(l100, l101, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz))
     *          + (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dbi)
     *          * (dmix(l110, l111, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz))
     *   dt/dtz = d(smoothStep(t).x)/dtz = 0
     *   therefore
     *   dz/dtz =
     *     (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/da)
     *     * ((dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dai) * (dmix(l000, l001, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)
     *        + (dmix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y)/dbi) * (dmix(l010, l011, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)))
     *   + (dmix(mix(mix(l000, l001, smoothStep(t).z), mix(l010, l011, smoothStep(t).z), smoothStep(t).y), mix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y), smoothStep(t).x)/db)
     *      * ((dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dai) * (dmix(l100, l101, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz))
     *         + (dmix(mix(l100, l101, smoothStep(t).z), mix(l110, l111, smoothStep(t).z), smoothStep(t).y)/dbi) * (dmix(l110, l111, smoothStep(t).z)/dtii) * (d(smoothStep(t).z)/dtz)))
     */

    /* Direct rewrite without the usage of temporary variables:
    return vec3(
        // dz/dtx
        mixDerivativeT(
            mix(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ),
            mix(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ),
            tMapped.x
        ) * tMappedDerivative.x,

        // dz/dty
        mixDerivativeA(
            mix(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ),
            mix(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ),
            tMapped.x
        ) * mixDerivativeT(
            mix(l000Hash, l001Hash, tMapped.z),
            mix(l010Hash, l011Hash, tMapped.z),
            tMapped.y
        ) * tMappedDerivative.y
        + mixDerivativeB(
            mix(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ),
            mix(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ),
            tMapped.x
        ) * mixDerivativeT(
            mix(l100Hash, l101Hash, tMapped.z),
            mix(l110Hash, l111Hash, tMapped.z),
            tMapped.y
        ) * tMappedDerivative.y,

        // dz/dtz
        mixDerivativeA(
            mix(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ),
            mix(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ),
            tMapped.x
        ) * (
            mixDerivativeA(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ) * mixDerivativeT(l000Hash, l001Hash, tMapped.z) * tMappedDerivative.z
            + mixDerivativeB(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ) * mixDerivativeT(l010Hash, l011Hash, tMapped.z) * tMappedDerivative.z
        ) + mixDerivativeB(
            mix(
                mix(l000Hash, l001Hash, tMapped.z),
                mix(l010Hash, l011Hash, tMapped.z),
                tMapped.y
            ),
            mix(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ),
            tMapped.x
        ) * (
            mixDerivativeA(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ) * mixDerivativeT(l100Hash, l101Hash, tMapped.z) * tMappedDerivative.z
            + mixDerivativeB(
                mix(l100Hash, l101Hash, tMapped.z),
                mix(l110Hash, l111Hash, tMapped.z),
                tMapped.y
            ) * mixDerivativeT(l110Hash, l111Hash, tMapped.z) * tMappedDerivative.z
        )
    */

    float mixl00XHash = mix(l000Hash, l001Hash, tMapped.z);
    float mixl01XHash = mix(l010Hash, l011Hash, tMapped.z);
    float mixl10XHash = mix(l100Hash, l101Hash, tMapped.z);
    float mixl11XHash = mix(l110Hash, l111Hash, tMapped.z);
    float mixl0XXHash = mix(mixl00XHash, mixl01XHash, tMapped.y);
    float mixl1XXHash = mix(mixl10XHash, mixl11XHash, tMapped.y);
    float mixDerivativeAlXXXHash = mixDerivativeA(mixl0XXHash, mixl1XXHash, tMapped.x);
    float mixDerivativeBlXXXHash = mixDerivativeB(mixl0XXHash, mixl1XXHash, tMapped.x);

    return vec3(
        // dz/dtx
        mixDerivativeT(
            mixl0XXHash,
            mixl1XXHash,
            tMapped.x
        ) * tMappedDerivative.x,

        // dz/dty
        (
            mixDerivativeAlXXXHash * mixDerivativeT(
                mixl00XHash,
                mixl01XHash,
                tMapped.y
            ) + mixDerivativeBlXXXHash * mixDerivativeT(
                mixl10XHash,
                mixl11XHash,
                tMapped.y
            )
        ) * tMappedDerivative.y,

        // dz/dtz
        (
            mixDerivativeAlXXXHash * (
                mixDerivativeA(
                    mixl00XHash,
                    mixl01XHash,
                    tMapped.y
                ) * mixDerivativeT(l000Hash, l001Hash, tMapped.z)
                + mixDerivativeB(
                    mixl00XHash,
                    mixl01XHash,
                    tMapped.y
                ) * mixDerivativeT(l010Hash, l011Hash, tMapped.z)
            ) + mixDerivativeBlXXXHash * (
                mixDerivativeA(
                    mixl10XHash,
                    mixl11XHash,
                    tMapped.y
                ) * mixDerivativeT(l100Hash, l101Hash, tMapped.z)
                + mixDerivativeB(
                    mixl10XHash,
                    mixl11XHash,
                    tMapped.y
                ) * mixDerivativeT(l110Hash, l111Hash, tMapped.z)
            )
        ) * tMappedDerivative.z
    );

    /* return vec2( */
    /*     mixDerivativeT( */
    /*         mix( */
    /*             mix(l000Hash, l001Hash, tMapped.z), */
    /*             mix(l010Hash, l011Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         mix( */
    /*             mix(l100Hash, l101Hash, tMapped.z), */
    /*             mix(l110Hash, l111Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         tMapped.x */
    /*     ) * tMappedDerivative.x, */

    /*     mixDerivativeA( */
    /*         mix( */
    /*             mix(l000Hash, l001Hash, tMapped.z), */
    /*             mix(l010Hash, l011Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         mix( */
    /*             mix(l100Hash, l101Hash, tMapped.z), */
    /*             mix(l110Hash, l111Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         tMapped.x */
    /*     ) * mixDerivativeT(l00Hash, l01Hash, tMapped.y) * tMappedDerivative.y */
    /*     + mixDerivativeB( */
    /*         mix( */
    /*             mix(l000Hash, l001Hash, tMapped.z), */
    /*             mix(l010Hash, l011Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         mix( */
    /*             mix(l100Hash, l101Hash, tMapped.z), */
    /*             mix(l110Hash, l111Hash, tMapped.z), */
    /*             tMapped.y */
    /*         ), */
    /*         tMapped.x */
    /*     ) * mixDerivativeT(l10Hash, l11Hash, tMapped.y) * tMappedDerivative.y */
    /* ); */

    /* return mix( */
    /*     mix( */
    /*         mix(l000Hash, l001Hash, tMapped.z), */
    /*         mix(l010Hash, l011Hash, tMapped.z), */
    /*         tMapped.y */
    /*     ), */
    /*     mix( */
    /*         mix(l100Hash, l101Hash, tMapped.z), */
    /*         mix(l110Hash, l111Hash, tMapped.z), */
    /*         tMapped.y */
    /*     ), */
    /*     tMapped.x */
    /* ); */
}

float valueNoiseUnsigned4(vec4 v) {
    // Lattice point integer coordinates
    ivec4 l0000 = ivec4(floor(v));
    ivec4 l0001 = l0000 + ivec4(0, 0, 0, 1);
    ivec4 l0010 = l0000 + ivec4(0, 0, 1, 0);
    ivec4 l0011 = l0000 + ivec4(0, 0, 1, 1);
    ivec4 l0100 = l0000 + ivec4(0, 1, 0, 0);
    ivec4 l0101 = l0000 + ivec4(0, 1, 0, 1);
    ivec4 l0110 = l0000 + ivec4(0, 1, 1, 0);
    ivec4 l0111 = l0000 + ivec4(0, 1, 1, 1);
    ivec4 l1000 = l0000 + ivec4(1, 0, 0, 0);
    ivec4 l1001 = l0000 + ivec4(1, 0, 0, 1);
    ivec4 l1010 = l0000 + ivec4(1, 0, 1, 0);
    ivec4 l1011 = l0000 + ivec4(1, 0, 1, 1);
    ivec4 l1100 = l0000 + ivec4(1, 1, 0, 0);
    ivec4 l1101 = l0000 + ivec4(1, 1, 0, 1);
    ivec4 l1110 = l0000 + ivec4(1, 1, 1, 0);
    ivec4 l1111 = l0000 + ivec4(1, 1, 1, 1);
    // Deterministically random values for lattice points
    float l0000Hash = ihash4AsUnsignedFloat(l0000);
    float l0001Hash = ihash4AsUnsignedFloat(l0001);
    float l0010Hash = ihash4AsUnsignedFloat(l0010);
    float l0011Hash = ihash4AsUnsignedFloat(l0011);
    float l0100Hash = ihash4AsUnsignedFloat(l0100);
    float l0101Hash = ihash4AsUnsignedFloat(l0101);
    float l0110Hash = ihash4AsUnsignedFloat(l0110);
    float l0111Hash = ihash4AsUnsignedFloat(l0111);
    float l1000Hash = ihash4AsUnsignedFloat(l1000);
    float l1001Hash = ihash4AsUnsignedFloat(l1001);
    float l1010Hash = ihash4AsUnsignedFloat(l1010);
    float l1011Hash = ihash4AsUnsignedFloat(l1011);
    float l1100Hash = ihash4AsUnsignedFloat(l1100);
    float l1101Hash = ihash4AsUnsignedFloat(l1101);
    float l1110Hash = ihash4AsUnsignedFloat(l1110);
    float l1111Hash = ihash4AsUnsignedFloat(l1111);
    // Interpolation
    vec4 t = fract(v);
    vec4 tMapped = noiseInterpolation(t);

    return mix(
        mix(
            mix(
                mix(l0000Hash, l0001Hash, tMapped.w),
                mix(l0010Hash, l0011Hash, tMapped.w),
                tMapped.z
            ),
            mix(
                mix(l0100Hash, l0101Hash, tMapped.w),
                mix(l0110Hash, l0111Hash, tMapped.w),
                tMapped.z
            ),
            tMapped.y
        ),
        mix(
            mix(
                mix(l1000Hash, l1001Hash, tMapped.w),
                mix(l1010Hash, l1011Hash, tMapped.w),
                tMapped.z
            ),
            mix(
                mix(l1100Hash, l1101Hash, tMapped.w),
                mix(l1110Hash, l1111Hash, tMapped.w),
                tMapped.z
            ),
            tMapped.y
        ),
        tMapped.x
    );
}

vec4 valueNoiseUnsignedGradient4(vec4 v) {
    // Lattice point integer coordinates
    ivec4 l0000 = ivec4(floor(v));
    ivec4 l0001 = l0000 + ivec4(0, 0, 0, 1);
    ivec4 l0010 = l0000 + ivec4(0, 0, 1, 0);
    ivec4 l0011 = l0000 + ivec4(0, 0, 1, 1);
    ivec4 l0100 = l0000 + ivec4(0, 1, 0, 0);
    ivec4 l0101 = l0000 + ivec4(0, 1, 0, 1);
    ivec4 l0110 = l0000 + ivec4(0, 1, 1, 0);
    ivec4 l0111 = l0000 + ivec4(0, 1, 1, 1);
    ivec4 l1000 = l0000 + ivec4(1, 0, 0, 0);
    ivec4 l1001 = l0000 + ivec4(1, 0, 0, 1);
    ivec4 l1010 = l0000 + ivec4(1, 0, 1, 0);
    ivec4 l1011 = l0000 + ivec4(1, 0, 1, 1);
    ivec4 l1100 = l0000 + ivec4(1, 1, 0, 0);
    ivec4 l1101 = l0000 + ivec4(1, 1, 0, 1);
    ivec4 l1110 = l0000 + ivec4(1, 1, 1, 0);
    ivec4 l1111 = l0000 + ivec4(1, 1, 1, 1);
    // Deterministically random values for lattice points
    float l0000Hash = ihash4AsUnsignedFloat(l0000);
    float l0001Hash = ihash4AsUnsignedFloat(l0001);
    float l0010Hash = ihash4AsUnsignedFloat(l0010);
    float l0011Hash = ihash4AsUnsignedFloat(l0011);
    float l0100Hash = ihash4AsUnsignedFloat(l0100);
    float l0101Hash = ihash4AsUnsignedFloat(l0101);
    float l0110Hash = ihash4AsUnsignedFloat(l0110);
    float l0111Hash = ihash4AsUnsignedFloat(l0111);
    float l1000Hash = ihash4AsUnsignedFloat(l1000);
    float l1001Hash = ihash4AsUnsignedFloat(l1001);
    float l1010Hash = ihash4AsUnsignedFloat(l1010);
    float l1011Hash = ihash4AsUnsignedFloat(l1011);
    float l1100Hash = ihash4AsUnsignedFloat(l1100);
    float l1101Hash = ihash4AsUnsignedFloat(l1101);
    float l1110Hash = ihash4AsUnsignedFloat(l1110);
    float l1111Hash = ihash4AsUnsignedFloat(l1111);
    // Interpolation
    vec4 t = fract(v);
    vec4 tMapped = noiseInterpolation(t);
    vec4 tMappedDerivative = noiseInterpolationGradient(t);

    float mixl000XHash = mix(l0000Hash, l0001Hash, tMapped.w);
    float mixl001XHash = mix(l0010Hash, l0011Hash, tMapped.w);
    float mixl010XHash = mix(l0100Hash, l0101Hash, tMapped.w);
    float mixl011XHash = mix(l0110Hash, l0111Hash, tMapped.w);
    float mixl100XHash = mix(l1000Hash, l1001Hash, tMapped.w);
    float mixl101XHash = mix(l1010Hash, l1011Hash, tMapped.w);
    float mixl110XHash = mix(l1100Hash, l1101Hash, tMapped.w);
    float mixl111XHash = mix(l1110Hash, l1111Hash, tMapped.w);
    float mixl00XXHash = mix(mixl000XHash, mixl001XHash, tMapped.z);
    float mixl01XXHash = mix(mixl010XHash, mixl011XHash, tMapped.z);
    float mixl10XXHash = mix(mixl100XHash, mixl101XHash, tMapped.z);
    float mixl11XXHash = mix(mixl110XHash, mixl111XHash, tMapped.z);
    float mixl0XXXHash = mix(mixl00XXHash, mixl01XXHash, tMapped.y);
    float mixl1XXXHash = mix(mixl10XXHash, mixl11XXHash, tMapped.y);
    float mixDerivativeAl0XXXHash = mixDerivativeA(mixl00XXHash, mixl01XXHash, tMapped.y);
    float mixDerivativeAl1XXXHash = mixDerivativeA(mixl10XXHash, mixl11XXHash, tMapped.y);
    float mixDerivativeBl0XXXHash = mixDerivativeB(mixl00XXHash, mixl01XXHash, tMapped.y);
    float mixDerivativeBl1XXXHash = mixDerivativeB(mixl10XXHash, mixl11XXHash, tMapped.y);
    float mixDerivativeAlXXXXHash = mixDerivativeA(mixl0XXXHash, mixl1XXXHash, tMapped.x);
    float mixDerivativeBlXXXXHash = mixDerivativeB(mixl0XXXHash, mixl1XXXHash, tMapped.x);

    // Pattern-extended to 4th dimension, without proof

    return vec4(
        // dz/dtx
        mixDerivativeT(
            mixl0XXXHash,
            mixl1XXXHash,
            tMapped.x
        ) * tMappedDerivative.x,

        // dz/dty
        (
            mixDerivativeAlXXXXHash
            * mixDerivativeT(mixl00XXHash, mixl01XXHash, tMapped.y)
            + mixDerivativeBlXXXXHash
            * mixDerivativeT(mixl10XXHash, mixl11XXHash, tMapped.y)
        ) * tMappedDerivative.y,

        // dz/dtz
        (
            mixDerivativeAlXXXXHash * (
                  mixDerivativeAl0XXXHash
                * mixDerivativeT(mixl000XHash, mixl001XHash, tMapped.z)
                + mixDerivativeBl0XXXHash
                * mixDerivativeT(mixl010XHash, mixl011XHash, tMapped.z)
            ) + mixDerivativeBlXXXXHash * (
                  mixDerivativeAl1XXXHash
                * mixDerivativeT(mixl100XHash, mixl101XHash, tMapped.z)
                + mixDerivativeBl1XXXHash
                * mixDerivativeT(mixl110XHash, mixl111XHash, tMapped.z)
            )
        ) * tMappedDerivative.z,

        // dz/dtw
        (
            mixDerivativeAlXXXXHash * (
                mixDerivativeAl0XXXHash * (
                      mixDerivativeA(mixl000XHash, mixl001XHash, tMapped.z)
                    * mixDerivativeT(   l0000Hash,    l0001Hash, tMapped.w)
                    + mixDerivativeB(mixl000XHash, mixl001XHash, tMapped.z)
                    * mixDerivativeT(   l0010Hash,    l0011Hash, tMapped.w)
                )
                + mixDerivativeBl0XXXHash * (
                      mixDerivativeA(mixl010XHash, mixl011XHash, tMapped.z)
                    * mixDerivativeT(   l0100Hash,    l0101Hash, tMapped.w)
                    + mixDerivativeB(mixl010XHash, mixl011XHash, tMapped.z)
                    * mixDerivativeT(   l0110Hash,    l0111Hash, tMapped.w)
                )
            ) + mixDerivativeBlXXXXHash * (
                mixDerivativeAl1XXXHash * (
                      mixDerivativeA(mixl100XHash, mixl101XHash, tMapped.z)
                    * mixDerivativeT(   l1000Hash,    l1001Hash, tMapped.w)
                    + mixDerivativeB(mixl100XHash, mixl101XHash, tMapped.z)
                    * mixDerivativeT(   l1010Hash,    l1011Hash, tMapped.w)
                )
                + mixDerivativeBl1XXXHash * (
                      mixDerivativeA(mixl110XHash, mixl111XHash, tMapped.z)
                    * mixDerivativeT(   l1100Hash,    l1101Hash, tMapped.w)
                    + mixDerivativeB(mixl110XHash, mixl111XHash, tMapped.z)
                    * mixDerivativeT(   l1110Hash,    l1111Hash, tMapped.w)
                )
            )
        ) * tMappedDerivative.w
    );
}

float valueNoiseSigned1(float x) { return valueNoiseUnsigned1(x) * 2.0 - 1.0; }
float valueNoiseSigned2(vec2 v)  { return valueNoiseUnsigned2(v) * 2.0 - 1.0; }
float valueNoiseSigned3(vec3 v)  { return valueNoiseUnsigned3(v) * 2.0 - 1.0; }
float valueNoiseSigned4(vec4 v)  { return valueNoiseUnsigned4(v) * 2.0 - 1.0; }
float valueNoiseSignedGradient1(float x) { return valueNoiseUnsignedGradient1(x) * 2.0; }
vec2  valueNoiseSignedGradient2(vec2 v)  { return valueNoiseUnsignedGradient2(v) * 2.0; }
vec3  valueNoiseSignedGradient3(vec3 v)  { return valueNoiseUnsignedGradient3(v) * 2.0; }
vec4  valueNoiseSignedGradient4(vec4 v)  { return valueNoiseUnsignedGradient4(v) * 2.0; }

float perlinNoiseSigned1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnitVec1(l0);
    float l1Hash = ihash1AsUnitVec1(l1);
    // Vectors to X from each lattice point
    float v0 = fract(x);
    float v1 = float(v0 - 1);
    // Interpolation
    float t = v0;
    float tMapped = noiseInterpolation(t);

    return mix(dot(l0Hash, v0), dot(l1Hash, v1), tMapped);
}

float perlinNoiseSignedGradient1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnitVec1(l0);
    float l1Hash = ihash1AsUnitVec1(l1);
    // Vectors to X from each lattice point
    float v0 = fract(x);
    float v1 = float(v0 - 1);
    // Interpolation
    float t = v0;
    float tMapped = noiseInterpolation(t);
    float tMappedDerivative = noiseInterpolationGradient(t);

    /*
        f(t) = mix(dot(l0Hash, v0), dot(l1Hash, v1), noiseInterpolation(t))
        df(t)/dt = dmix(dot(l0Hash, v0), dot(l1Hash, v1), noiseInterpolation(t))/dt
            = (dmix(dot(l0Hash, v0), dot(l1Hash, v1), noiseInterpolation(t))/da) * (da/dt)
            + (dmix(dot(l0Hash, v0), dot(l1Hash, v1), noiseInterpolation(t))/db) * (db/dt)
            + (dmix(dot(l0Hash, v0), dot(l1Hash, v1), noiseInterpolation(t))/dt2) * (dt2/dt)
        where
          dot(x, y) = x * y
          ddot(x, y)/dx = d(x * y)/dx = y
          ddot(x, y)/dy = d(x * y)/dy = x
          da/dt  = ddot(l0Hash, v0)/dt = ddot(l0Hash,     t)/dt
                 = (ddot(l0Hash,     t)/dx) * (dx/dt) (= 0)
                 + (ddot(l0Hash,     t)/dy) * (dy/dt)
          where
              dx/dt = dl0Hash/dt = 0
              dy/dt = dt/dt = 1
          therefore
              da/dt = (ddot(l0Hash,     t)/dy) * 1 = l0Hash
          db/dt  = ddot(l1Hash, v1)/dt = ddot(l1Hash, 1 - t)/dt
                 = (ddot(l1Hash, 1 - t)/dx) * (dx/dt) (= 0)
                 + (ddot(l1Hash, 1 - t)/dy) * (dy/dt)
          where
              dx/dt = dl1Hash/dt = 0
              dy/dt = d(1 - t)/dt = -1
          therefore
              db/dt = (ddot(l1Hash, 1 - t)/dy) * (-1) = -l1Hash
          dt2/dt = dnoiseInterpolation(t)/dt = noiseInterpolationGradient(t)
     */

    float dot0 = dot(l0Hash, v0);
    float dot1 = dot(l1Hash, v1);

    return mixDerivativeA(dot0, dot1, tMapped) * l0Hash
         + mixDerivativeB(dot0, dot1, tMapped) * l1Hash
         + mixDerivativeT(dot0, dot1, tMapped) * tMappedDerivative;
}

float perlinNoiseSigned2(vec2 v) {
    // Lattice point integer coordinates
    ivec2 l00 = ivec2(floor(v));
    ivec2 l01 = l00 + ivec2(0, 1);
    ivec2 l10 = l00 + ivec2(1, 0);
    ivec2 l11 = l00 + ivec2(1, 1);
    // Deterministically random values for lattice points
    vec2 l00Hash = ihash2AsUnitVec2(l00);
    vec2 l01Hash = ihash2AsUnitVec2(l01);
    vec2 l10Hash = ihash2AsUnitVec2(l10);
    vec2 l11Hash = ihash2AsUnitVec2(l11);
    // Vectors to V from each lattice point
    vec2 v00 = fract(v);
    vec2 v01 = vec2(v00.x    , v00.y - 1);
    vec2 v10 = vec2(v00.x - 1, v00.y    );
    vec2 v11 = vec2(v00.x - 1, v00.y - 1);
    // Interpolation
    vec2 t = v00;
    vec2 tMapped = noiseInterpolation(t);

    return mix(
        mix(
            dot(l00Hash, v00),
            dot(l01Hash, v01),
            tMapped.y
        ),
        mix(
            dot(l10Hash, v10),
            dot(l11Hash, v11),
            tMapped.y
        ),
        tMapped.x
    );
}

vec2 perlinNoiseSignedGradient2(vec2 v) {
    // Lattice point integer coordinates
    ivec2 l00 = ivec2(floor(v));
    ivec2 l01 = l00 + ivec2(0, 1);
    ivec2 l10 = l00 + ivec2(1, 0);
    ivec2 l11 = l00 + ivec2(1, 1);
    // Deterministically random values for lattice points
    vec2 l00Hash = ihash2AsUnitVec2(l00);
    vec2 l01Hash = ihash2AsUnitVec2(l01);
    vec2 l10Hash = ihash2AsUnitVec2(l10);
    vec2 l11Hash = ihash2AsUnitVec2(l11);
    // Vectors to V from each lattice point
    vec2 v00 = fract(v);
    vec2 v01 = vec2(v00.x    , v00.y - 1);
    vec2 v10 = vec2(v00.x - 1, v00.y    );
    vec2 v11 = vec2(v00.x - 1, v00.y - 1);
    // Interpolation
    vec2 t = v00;
    vec2 tMapped = noiseInterpolation(t);
    vec2 tMappedDerivative = noiseInterpolationGradient(t);

    /*
    f(t) = mix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))
    df(t)/dtx = dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dtx
    df(t)/dtx = (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/da) * (da/dtx)
              + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/db) * (db/dtx)
              + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dt2) * (dt2/dtx)
    where
        da/dtx  = dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dtx
                = (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) * (dai/dtx)
                + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * (dbi/dtx)
                + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dti) * (dti/dtx)
        where
            dot(i, j) = ix * jx + iy * jy
            ddot(i, j)/dix = d(ix * jx + iy * jy)/dix = jx
            dai/dtx = ddot(l00Hash, v00)/dtx = ddot(l00Hash, tx)/dtx = l00Hash.x
            dbi/dtx = ddot(l01Hash, v01)/dtx = ddot(l01Hash, tx)/dtx = l01Hash.x
            dti/dtx = dinterp(ty)/dtx = 0
        therefore
            da/dtx = (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) * l00Hash.x
                   + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * l01Hash.x
        db/dtx = (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dai) * (-l10Hash.x)
               + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dbi) * (-l11Hash.x)
        dt2/dtx = dinterp(tx)/dtx
    therefore
        df(t)/dtx = (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/da) * (
              (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) * l00Hash.x
            + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * l01Hash.x
        ) + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/db) * (
              (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dai) * (-l10Hash.x)
            + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dbi) * (-l11Hash.x)
        ) + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dt2)
             * (dinterp(tx)/dtx)
        //df(t)/dtx = (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dt2) * dinterp(tx)/dtx

    df(t)/dty = dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dty
    df(t)/dty = (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/da) * (da/dty)
              + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/db) * (db/dty)
              + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/dt2) * (dt2/dty)
    where
        da/dty  = dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dty =
                = (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) * (dai/dty)
                + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * (dbi/dty)
                + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dti) * (dti/dty)
        where
            dai/dty = ddot(l00Hash, v00)/dty =  l01Hash.y
            dbi/dty = ddot(l01Hash, v01)/dty = -l01Hash.y
            dti/dty = dinterp(ty)/dty
        therefore
            da/dty = (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) *   l00Hash.y
                   + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * (-l01Hash.y)
                   + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dti) * (dinterp(ty)/dty)
        db/dty = (analogously) = (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dai) *   l10Hash.y
                               + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dbi) * (-l11Hash.y)
                               + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dti) * (dinterp(ty)/dty)
        dt2/dty = dinterp(tx)/dty = 0
    therefore
        df(t)/dty = (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/da) * (
                        (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dai) *   l00Hash.y
                        + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dbi) * (-l01Hash.y)
                        + (dmix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty))/dti) * (dinterp(ty)/dty)
                    ) + (dmix(mix(dot(l00Hash, v00), dot(l01Hash, v01), interp(ty)), mix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty)), interp(tx))/db) * (
                        (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dai) *   l10Hash.y
                        + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dbi) * (-l11Hash.y)
                        + (dmix(dot(l10Hash, v10), dot(l11Hash, v11), interp(ty))/dti) * (dinterp(ty)/dty)
                    )
    */

    float dot00 = dot(l00Hash, v00);
    float dot01 = dot(l01Hash, v01);
    float dot10 = dot(l10Hash, v10);
    float dot11 = dot(l11Hash, v11);
    float mixdot0X = mix(dot00, dot01, tMapped.y);
    float mixdot1X = mix(dot10, dot11, tMapped.y);
    float mixDerivativeAdot0X = mixDerivativeA(dot00, dot01, tMapped.y);
    float mixDerivativeBdot0X = mixDerivativeB(dot00, dot01, tMapped.y);
    float mixDerivativeAdot1X = mixDerivativeA(dot10, dot11, tMapped.y);
    float mixDerivativeBdot1X = mixDerivativeB(dot10, dot11, tMapped.y);
    float mixDerivativeAdotXX = mixDerivativeA(mixdot0X, mixdot1X, tMapped.x);
    float mixDerivativeBdotXX = mixDerivativeB(mixdot0X, mixdot1X, tMapped.x);

    return vec2(
        mixDerivativeAdotXX * (
              mixDerivativeAdot0X * l00Hash.x
            + mixDerivativeBdot0X * l01Hash.x
        ) + mixDerivativeBdotXX * (
              mixDerivativeAdot1X * l10Hash.x
            + mixDerivativeBdot1X * l11Hash.x
        ) + mixDerivativeT(mixdot0X, mixdot1X, tMapped.x) * tMappedDerivative.x,

        mixDerivativeAdotXX * (
             mixDerivativeAdot0X * l00Hash.y
           + mixDerivativeBdot0X * l01Hash.y
           + mixDerivativeT(dot00, dot01, tMapped.y) * tMappedDerivative.y
        ) + mixDerivativeBdotXX * (
             mixDerivativeAdot1X * l10Hash.y
           + mixDerivativeBdot1X * l11Hash.y
           + mixDerivativeT(dot10, dot11, tMapped.y) * tMappedDerivative.y
        )
    );
}

float perlinNoiseSigned3(vec3 v) {
    // Lattice point integer coordinates
    ivec3 l000 = ivec3(floor(v));
    ivec3 l001 = l000 + ivec3(0, 0, 1);
    ivec3 l010 = l000 + ivec3(0, 1, 0);
    ivec3 l011 = l000 + ivec3(0, 1, 1);
    ivec3 l100 = l000 + ivec3(1, 0, 0);
    ivec3 l101 = l000 + ivec3(1, 0, 1);
    ivec3 l110 = l000 + ivec3(1, 1, 0);
    ivec3 l111 = l000 + ivec3(1, 1, 1);
    // Deterministically random values for lattice points
    vec3 l000Hash = ihash3AsUnitVec3(l000);
    vec3 l001Hash = ihash3AsUnitVec3(l001);
    vec3 l010Hash = ihash3AsUnitVec3(l010);
    vec3 l011Hash = ihash3AsUnitVec3(l011);
    vec3 l100Hash = ihash3AsUnitVec3(l100);
    vec3 l101Hash = ihash3AsUnitVec3(l101);
    vec3 l110Hash = ihash3AsUnitVec3(l110);
    vec3 l111Hash = ihash3AsUnitVec3(l111);
    // Vectors to V from each lattice point
    vec3 v000 = fract(v);
    vec3 v001 = vec3(v000.x    , v000.y    , v000.z - 1);
    vec3 v010 = vec3(v000.x    , v000.y - 1, v000.z    );
    vec3 v011 = vec3(v000.x    , v000.y - 1, v000.z - 1);
    vec3 v100 = vec3(v000.x - 1, v000.y    , v000.z    );
    vec3 v101 = vec3(v000.x - 1, v000.y    , v000.z - 1);
    vec3 v110 = vec3(v000.x - 1, v000.y - 1, v000.z    );
    vec3 v111 = vec3(v000.x - 1, v000.y - 1, v000.z - 1);
    // Interpolation
    vec3 t = v000;
    vec3 tMapped = noiseInterpolation(t);

    return mix(
        mix(
            mix(
                dot(l000Hash, v000),
                dot(l001Hash, v001),
                tMapped.z
            ),
            mix(
                dot(l010Hash, v010),
                dot(l011Hash, v011),
                tMapped.z
            ),
            tMapped.y
        ),
        mix(
            mix(
                dot(l100Hash, v100),
                dot(l101Hash, v101),
                tMapped.z
            ),
            mix(
                dot(l110Hash, v110),
                dot(l111Hash, v111),
                tMapped.z
            ),
            tMapped.y
        ),
        tMapped.x
    );
}

vec3 perlinNoiseSignedGradient3(vec3 v) {
    // Lattice point integer coordinates
    ivec3 l000 = ivec3(floor(v));
    ivec3 l001 = l000 + ivec3(0, 0, 1);
    ivec3 l010 = l000 + ivec3(0, 1, 0);
    ivec3 l011 = l000 + ivec3(0, 1, 1);
    ivec3 l100 = l000 + ivec3(1, 0, 0);
    ivec3 l101 = l000 + ivec3(1, 0, 1);
    ivec3 l110 = l000 + ivec3(1, 1, 0);
    ivec3 l111 = l000 + ivec3(1, 1, 1);
    // Deterministically random values for lattice points
    vec3 l000Hash = ihash3AsUnitVec3(l000);
    vec3 l001Hash = ihash3AsUnitVec3(l001);
    vec3 l010Hash = ihash3AsUnitVec3(l010);
    vec3 l011Hash = ihash3AsUnitVec3(l011);
    vec3 l100Hash = ihash3AsUnitVec3(l100);
    vec3 l101Hash = ihash3AsUnitVec3(l101);
    vec3 l110Hash = ihash3AsUnitVec3(l110);
    vec3 l111Hash = ihash3AsUnitVec3(l111);
    // Vectors to V from each lattice point
    vec3 v000 = fract(v);
    vec3 v001 = vec3(v000.x    , v000.y    , v000.z - 1);
    vec3 v010 = vec3(v000.x    , v000.y - 1, v000.z    );
    vec3 v011 = vec3(v000.x    , v000.y - 1, v000.z - 1);
    vec3 v100 = vec3(v000.x - 1, v000.y    , v000.z    );
    vec3 v101 = vec3(v000.x - 1, v000.y    , v000.z - 1);
    vec3 v110 = vec3(v000.x - 1, v000.y - 1, v000.z    );
    vec3 v111 = vec3(v000.x - 1, v000.y - 1, v000.z - 1);
    // Interpolation
    vec3 t = v000;
    vec3 tMapped = noiseInterpolation(t);
    vec3 tMappedDerivative = noiseInterpolationGradient(t);

    // Pattern-expanded from the previous dimension, no proof given

    float dot000 = dot(l000Hash, v000);
    float dot001 = dot(l001Hash, v001);
    float dot010 = dot(l010Hash, v010);
    float dot011 = dot(l011Hash, v011);
    float dot100 = dot(l100Hash, v100);
    float dot101 = dot(l101Hash, v101);
    float dot110 = dot(l110Hash, v110);
    float dot111 = dot(l111Hash, v111);
    float mixdot00X = mix(dot000, dot001, tMapped.z);
    float mixdot01X = mix(dot010, dot011, tMapped.z);
    float mixdot10X = mix(dot100, dot101, tMapped.z);
    float mixdot11X = mix(dot110, dot111, tMapped.z);
    float mixdot0XX = mix(mixdot00X, mixdot01X, tMapped.y);
    float mixdot1XX = mix(mixdot10X, mixdot11X, tMapped.y);
    float mixDerivativeAdot00X = mixDerivativeA(dot000, dot001, tMapped.z);
    float mixDerivativeBdot00X = mixDerivativeB(dot000, dot001, tMapped.z);
    float mixDerivativeAdot01X = mixDerivativeA(dot010, dot011, tMapped.z);
    float mixDerivativeBdot01X = mixDerivativeB(dot010, dot011, tMapped.z);
    float mixDerivativeAdot10X = mixDerivativeA(dot100, dot101, tMapped.z);
    float mixDerivativeBdot10X = mixDerivativeB(dot100, dot101, tMapped.z);
    float mixDerivativeAdot11X = mixDerivativeA(dot110, dot111, tMapped.z);
    float mixDerivativeBdot11X = mixDerivativeB(dot110, dot111, tMapped.z);
    float mixDerivativeAdot0XX = mixDerivativeA(mixdot00X, mixdot01X, tMapped.y);
    float mixDerivativeBdot0XX = mixDerivativeB(mixdot00X, mixdot01X, tMapped.y);
    float mixDerivativeAdot1XX = mixDerivativeA(mixdot10X, mixdot11X, tMapped.y);
    float mixDerivativeBdot1XX = mixDerivativeB(mixdot10X, mixdot11X, tMapped.y);
    float mixDerivativeAdotXXX = mixDerivativeA(mixdot0XX, mixdot1XX, tMapped.x);
    float mixDerivativeBdotXXX = mixDerivativeB(mixdot0XX, mixdot1XX, tMapped.x);

    return vec3(
        mixDerivativeAdotXXX * (
            mixDerivativeAdot0XX * (
                mixDerivativeAdot00X * l000Hash.x
              + mixDerivativeBdot00X * l001Hash.x
            )
            + mixDerivativeBdot0XX * (
                mixDerivativeAdot01X * l010Hash.x
              + mixDerivativeBdot01X * l011Hash.x
            )
        ) + mixDerivativeBdotXXX * (
            mixDerivativeAdot1XX * (
                mixDerivativeAdot10X * l100Hash.x
              + mixDerivativeBdot10X * l101Hash.x
            )
            + mixDerivativeBdot1XX * (
                mixDerivativeAdot11X * l110Hash.x
              + mixDerivativeBdot11X * l111Hash.x
            )
        ) + mixDerivativeT(mixdot0XX, mixdot1XX, tMapped.x) * tMappedDerivative.x,

        mixDerivativeAdotXXX * (
            mixDerivativeAdot0XX * (
                mixDerivativeAdot00X * l000Hash.y
              + mixDerivativeBdot00X * l001Hash.y
            )
            + mixDerivativeBdot0XX * (
                mixDerivativeAdot01X * l010Hash.y
              + mixDerivativeBdot01X * l011Hash.y
            )
            + mixDerivativeT(mixdot00X, mixdot01X, tMapped.y) * tMappedDerivative.y
        ) + mixDerivativeBdotXXX * (
            mixDerivativeAdot1XX * (
                mixDerivativeAdot10X * l100Hash.y
              + mixDerivativeBdot10X * l101Hash.y
            )
            + mixDerivativeBdot1XX * (
                mixDerivativeAdot11X * l110Hash.y
              + mixDerivativeBdot11X * l111Hash.y
            )
            + mixDerivativeT(mixdot10X, mixdot11X, tMapped.y) * tMappedDerivative.y
        ),

        mixDerivativeAdotXXX * (
            mixDerivativeAdot0XX * (
                mixDerivativeAdot00X * l000Hash.z
              + mixDerivativeBdot00X * l001Hash.z
              + mixDerivativeT(dot000, dot001, tMapped.z) * tMappedDerivative.z
            )
            + mixDerivativeBdot0XX * (
                mixDerivativeAdot01X * l010Hash.z
              + mixDerivativeBdot01X * l011Hash.z
              + mixDerivativeT(dot010, dot011, tMapped.z) * tMappedDerivative.z
            )
        ) + mixDerivativeBdotXXX * (
            mixDerivativeAdot1XX * (
                mixDerivativeAdot10X * l100Hash.z
              + mixDerivativeBdot10X * l101Hash.z
              + mixDerivativeT(dot100, dot101, tMapped.z) * tMappedDerivative.z
            )
            + mixDerivativeBdot1XX * (
                mixDerivativeAdot11X * l110Hash.z
              + mixDerivativeBdot11X * l111Hash.z
              + mixDerivativeT(dot110, dot111, tMapped.z) * tMappedDerivative.z
            )
        )
    );
}

float perlinNoiseSigned4(vec4 v) {
    // Lattice point integer coordinates
    ivec4 l0000 = ivec4(floor(v));
    ivec4 l0001 = l0000 + ivec4(0, 0, 0, 1);
    ivec4 l0010 = l0000 + ivec4(0, 0, 1, 0);
    ivec4 l0011 = l0000 + ivec4(0, 0, 1, 1);
    ivec4 l0100 = l0000 + ivec4(0, 1, 0, 0);
    ivec4 l0101 = l0000 + ivec4(0, 1, 0, 1);
    ivec4 l0110 = l0000 + ivec4(0, 1, 1, 0);
    ivec4 l0111 = l0000 + ivec4(0, 1, 1, 1);
    ivec4 l1000 = l0000 + ivec4(1, 0, 0, 0);
    ivec4 l1001 = l0000 + ivec4(1, 0, 0, 1);
    ivec4 l1010 = l0000 + ivec4(1, 0, 1, 0);
    ivec4 l1011 = l0000 + ivec4(1, 0, 1, 1);
    ivec4 l1100 = l0000 + ivec4(1, 1, 0, 0);
    ivec4 l1101 = l0000 + ivec4(1, 1, 0, 1);
    ivec4 l1110 = l0000 + ivec4(1, 1, 1, 0);
    ivec4 l1111 = l0000 + ivec4(1, 1, 1, 1);
    // Deterministically random values for lattice points
    vec4 l0000Hash = ihash4AsUnitVec4(l0000);
    vec4 l0001Hash = ihash4AsUnitVec4(l0001);
    vec4 l0010Hash = ihash4AsUnitVec4(l0010);
    vec4 l0011Hash = ihash4AsUnitVec4(l0011);
    vec4 l0100Hash = ihash4AsUnitVec4(l0100);
    vec4 l0101Hash = ihash4AsUnitVec4(l0101);
    vec4 l0110Hash = ihash4AsUnitVec4(l0110);
    vec4 l0111Hash = ihash4AsUnitVec4(l0111);
    vec4 l1000Hash = ihash4AsUnitVec4(l1000);
    vec4 l1001Hash = ihash4AsUnitVec4(l1001);
    vec4 l1010Hash = ihash4AsUnitVec4(l1010);
    vec4 l1011Hash = ihash4AsUnitVec4(l1011);
    vec4 l1100Hash = ihash4AsUnitVec4(l1100);
    vec4 l1101Hash = ihash4AsUnitVec4(l1101);
    vec4 l1110Hash = ihash4AsUnitVec4(l1110);
    vec4 l1111Hash = ihash4AsUnitVec4(l1111);
    // Vectors to V from each lattice point
    vec4 v0000 = fract(v);
    vec4 v0001 = vec4(v0000.x    , v0000.y    , v0000.z    , v0000.w - 1);
    vec4 v0010 = vec4(v0000.x    , v0000.y    , v0000.z - 1, v0000.w    );
    vec4 v0011 = vec4(v0000.x    , v0000.y    , v0000.z - 1, v0000.w - 1);
    vec4 v0100 = vec4(v0000.x    , v0000.y - 1, v0000.z    , v0000.w    );
    vec4 v0101 = vec4(v0000.x    , v0000.y - 1, v0000.z    , v0000.w - 1);
    vec4 v0110 = vec4(v0000.x    , v0000.y - 1, v0000.z - 1, v0000.w    );
    vec4 v0111 = vec4(v0000.x    , v0000.y - 1, v0000.z - 1, v0000.w - 1);
    vec4 v1000 = vec4(v0000.x - 1, v0000.y    , v0000.z    , v0000.w    );
    vec4 v1001 = vec4(v0000.x - 1, v0000.y    , v0000.z    , v0000.w - 1);
    vec4 v1010 = vec4(v0000.x - 1, v0000.y    , v0000.z - 1, v0000.w    );
    vec4 v1011 = vec4(v0000.x - 1, v0000.y    , v0000.z - 1, v0000.w - 1);
    vec4 v1100 = vec4(v0000.x - 1, v0000.y - 1, v0000.z    , v0000.w    );
    vec4 v1101 = vec4(v0000.x - 1, v0000.y - 1, v0000.z    , v0000.w - 1);
    vec4 v1110 = vec4(v0000.x - 1, v0000.y - 1, v0000.z - 1, v0000.w    );
    vec4 v1111 = vec4(v0000.x - 1, v0000.y - 1, v0000.z - 1, v0000.w - 1);
    // Interpolation
    vec4 t = v0000;
    vec4 tMapped = noiseInterpolation(t);

    return mix(
        mix(
            mix(
                mix(
                    dot(l0000Hash, v0000),
                    dot(l0001Hash, v0001),
                    tMapped.w
                ),
                mix(
                    dot(l0010Hash, v0010),
                    dot(l0011Hash, v0011),
                    tMapped.w
                ),
                tMapped.z
            ),
            mix(
                mix(
                    dot(l0100Hash, v0100),
                    dot(l0101Hash, v0101),
                    tMapped.w
                ),
                mix(
                    dot(l0110Hash, v0110),
                    dot(l0111Hash, v0111),
                    tMapped.w
                ),
                tMapped.z
            ),
            tMapped.y
        ),
        mix(
            mix(
                mix(
                    dot(l1000Hash, v1000),
                    dot(l1001Hash, v1001),
                    tMapped.w
                ),
                mix(
                    dot(l1010Hash, v1010),
                    dot(l1011Hash, v1011),
                    tMapped.w
                ),
                tMapped.z
            ),
            mix(
                mix(
                    dot(l1100Hash, v1100),
                    dot(l1101Hash, v1101),
                    tMapped.w
                ),
                mix(
                    dot(l1110Hash, v1110),
                    dot(l1111Hash, v1111),
                    tMapped.w
                ),
                tMapped.z
            ),
            tMapped.y
        ),
        tMapped.x
    );
}

vec4 perlinNoiseSignedGradient4(vec4 v) {
    // Lattice point integer coordinates
    ivec4 l0000 = ivec4(floor(v));
    ivec4 l0001 = l0000 + ivec4(0, 0, 0, 1);
    ivec4 l0010 = l0000 + ivec4(0, 0, 1, 0);
    ivec4 l0011 = l0000 + ivec4(0, 0, 1, 1);
    ivec4 l0100 = l0000 + ivec4(0, 1, 0, 0);
    ivec4 l0101 = l0000 + ivec4(0, 1, 0, 1);
    ivec4 l0110 = l0000 + ivec4(0, 1, 1, 0);
    ivec4 l0111 = l0000 + ivec4(0, 1, 1, 1);
    ivec4 l1000 = l0000 + ivec4(1, 0, 0, 0);
    ivec4 l1001 = l0000 + ivec4(1, 0, 0, 1);
    ivec4 l1010 = l0000 + ivec4(1, 0, 1, 0);
    ivec4 l1011 = l0000 + ivec4(1, 0, 1, 1);
    ivec4 l1100 = l0000 + ivec4(1, 1, 0, 0);
    ivec4 l1101 = l0000 + ivec4(1, 1, 0, 1);
    ivec4 l1110 = l0000 + ivec4(1, 1, 1, 0);
    ivec4 l1111 = l0000 + ivec4(1, 1, 1, 1);
    // Deterministically random values for lattice points
    vec4 l0000Hash = ihash4AsUnitVec4(l0000);
    vec4 l0001Hash = ihash4AsUnitVec4(l0001);
    vec4 l0010Hash = ihash4AsUnitVec4(l0010);
    vec4 l0011Hash = ihash4AsUnitVec4(l0011);
    vec4 l0100Hash = ihash4AsUnitVec4(l0100);
    vec4 l0101Hash = ihash4AsUnitVec4(l0101);
    vec4 l0110Hash = ihash4AsUnitVec4(l0110);
    vec4 l0111Hash = ihash4AsUnitVec4(l0111);
    vec4 l1000Hash = ihash4AsUnitVec4(l1000);
    vec4 l1001Hash = ihash4AsUnitVec4(l1001);
    vec4 l1010Hash = ihash4AsUnitVec4(l1010);
    vec4 l1011Hash = ihash4AsUnitVec4(l1011);
    vec4 l1100Hash = ihash4AsUnitVec4(l1100);
    vec4 l1101Hash = ihash4AsUnitVec4(l1101);
    vec4 l1110Hash = ihash4AsUnitVec4(l1110);
    vec4 l1111Hash = ihash4AsUnitVec4(l1111);
    // Vectors to V from each lattice point
    vec4 v0000 = fract(v);
    vec4 v0001 = vec4(v0000.x    , v0000.y    , v0000.z    , v0000.w - 1);
    vec4 v0010 = vec4(v0000.x    , v0000.y    , v0000.z - 1, v0000.w    );
    vec4 v0011 = vec4(v0000.x    , v0000.y    , v0000.z - 1, v0000.w - 1);
    vec4 v0100 = vec4(v0000.x    , v0000.y - 1, v0000.z    , v0000.w    );
    vec4 v0101 = vec4(v0000.x    , v0000.y - 1, v0000.z    , v0000.w - 1);
    vec4 v0110 = vec4(v0000.x    , v0000.y - 1, v0000.z - 1, v0000.w    );
    vec4 v0111 = vec4(v0000.x    , v0000.y - 1, v0000.z - 1, v0000.w - 1);
    vec4 v1000 = vec4(v0000.x - 1, v0000.y    , v0000.z    , v0000.w    );
    vec4 v1001 = vec4(v0000.x - 1, v0000.y    , v0000.z    , v0000.w - 1);
    vec4 v1010 = vec4(v0000.x - 1, v0000.y    , v0000.z - 1, v0000.w    );
    vec4 v1011 = vec4(v0000.x - 1, v0000.y    , v0000.z - 1, v0000.w - 1);
    vec4 v1100 = vec4(v0000.x - 1, v0000.y - 1, v0000.z    , v0000.w    );
    vec4 v1101 = vec4(v0000.x - 1, v0000.y - 1, v0000.z    , v0000.w - 1);
    vec4 v1110 = vec4(v0000.x - 1, v0000.y - 1, v0000.z - 1, v0000.w    );
    vec4 v1111 = vec4(v0000.x - 1, v0000.y - 1, v0000.z - 1, v0000.w - 1);
    // Interpolation
    vec4 t = v0000;
    vec4 tMapped = noiseInterpolation(t);
    vec4 tMappedDerivative = noiseInterpolationGradient(t);

    // Pattern-expanded from the previous dimension, no proof given

    float dot0000 = dot(l0000Hash, v0000);
    float dot0001 = dot(l0001Hash, v0001);
    float dot0010 = dot(l0010Hash, v0010);
    float dot0011 = dot(l0011Hash, v0011);
    float dot0100 = dot(l0100Hash, v0100);
    float dot0101 = dot(l0101Hash, v0101);
    float dot0110 = dot(l0110Hash, v0110);
    float dot0111 = dot(l0111Hash, v0111);
    float dot1000 = dot(l1000Hash, v1000);
    float dot1001 = dot(l1001Hash, v1001);
    float dot1010 = dot(l1010Hash, v1010);
    float dot1011 = dot(l1011Hash, v1011);
    float dot1100 = dot(l1100Hash, v1100);
    float dot1101 = dot(l1101Hash, v1101);
    float dot1110 = dot(l1110Hash, v1110);
    float dot1111 = dot(l1111Hash, v1111);
    float mixdot000X = mix(dot0000, dot0001, tMapped.w);
    float mixdot001X = mix(dot0010, dot0011, tMapped.w);
    float mixdot010X = mix(dot0100, dot0101, tMapped.w);
    float mixdot011X = mix(dot0110, dot0111, tMapped.w);
    float mixdot100X = mix(dot1000, dot1001, tMapped.w);
    float mixdot101X = mix(dot1010, dot1011, tMapped.w);
    float mixdot110X = mix(dot1100, dot1101, tMapped.w);
    float mixdot111X = mix(dot1110, dot1111, tMapped.w);
    float mixdot00XX = mix(mixdot000X, mixdot001X, tMapped.z);
    float mixdot01XX = mix(mixdot010X, mixdot011X, tMapped.z);
    float mixdot10XX = mix(mixdot100X, mixdot101X, tMapped.z);
    float mixdot11XX = mix(mixdot110X, mixdot111X, tMapped.z);
    float mixdot0XXX = mix(mixdot00XX, mixdot01XX, tMapped.y);
    float mixdot1XXX = mix(mixdot10XX, mixdot11XX, tMapped.y);
    float mixDerivativeAdot000X = mixDerivativeA(dot0000, dot0001, tMapped.w);
    float mixDerivativeBdot000X = mixDerivativeB(dot0000, dot0001, tMapped.w);
    float mixDerivativeAdot001X = mixDerivativeA(dot0010, dot0011, tMapped.w);
    float mixDerivativeBdot001X = mixDerivativeB(dot0010, dot0011, tMapped.w);
    float mixDerivativeAdot010X = mixDerivativeA(dot0100, dot0101, tMapped.w);
    float mixDerivativeBdot010X = mixDerivativeB(dot0100, dot0101, tMapped.w);
    float mixDerivativeAdot011X = mixDerivativeA(dot0110, dot0111, tMapped.w);
    float mixDerivativeBdot011X = mixDerivativeB(dot0110, dot0111, tMapped.w);
    float mixDerivativeAdot100X = mixDerivativeA(dot1000, dot1001, tMapped.w);
    float mixDerivativeBdot100X = mixDerivativeB(dot1000, dot1001, tMapped.w);
    float mixDerivativeAdot101X = mixDerivativeA(dot1010, dot1011, tMapped.w);
    float mixDerivativeBdot101X = mixDerivativeB(dot1010, dot1011, tMapped.w);
    float mixDerivativeAdot110X = mixDerivativeA(dot1100, dot1101, tMapped.w);
    float mixDerivativeBdot110X = mixDerivativeB(dot1100, dot1101, tMapped.w);
    float mixDerivativeAdot111X = mixDerivativeA(dot1110, dot1111, tMapped.w);
    float mixDerivativeBdot111X = mixDerivativeB(dot1110, dot1111, tMapped.w);
    float mixDerivativeAdot00XX = mixDerivativeA(mixdot000X, mixdot001X, tMapped.z);
    float mixDerivativeBdot00XX = mixDerivativeB(mixdot000X, mixdot001X, tMapped.z);
    float mixDerivativeAdot01XX = mixDerivativeA(mixdot010X, mixdot011X, tMapped.z);
    float mixDerivativeBdot01XX = mixDerivativeB(mixdot010X, mixdot011X, tMapped.z);
    float mixDerivativeAdot10XX = mixDerivativeA(mixdot100X, mixdot101X, tMapped.z);
    float mixDerivativeBdot10XX = mixDerivativeB(mixdot100X, mixdot101X, tMapped.z);
    float mixDerivativeAdot11XX = mixDerivativeA(mixdot110X, mixdot111X, tMapped.z);
    float mixDerivativeBdot11XX = mixDerivativeB(mixdot110X, mixdot111X, tMapped.z);
    float mixDerivativeAdot0XXX = mixDerivativeA(mixdot00XX, mixdot01XX, tMapped.y);
    float mixDerivativeBdot0XXX = mixDerivativeB(mixdot00XX, mixdot01XX, tMapped.y);
    float mixDerivativeAdot1XXX = mixDerivativeA(mixdot10XX, mixdot11XX, tMapped.y);
    float mixDerivativeBdot1XXX = mixDerivativeB(mixdot10XX, mixdot11XX, tMapped.y);
    float mixDerivativeAdotXXXX = mixDerivativeA(mixdot0XXX, mixdot1XXX, tMapped.x);
    float mixDerivativeBdotXXXX = mixDerivativeB(mixdot0XXX, mixdot1XXX, tMapped.x);

    return vec4(
        // dz/dtx
        mixDerivativeAdotXXXX * (
            mixDerivativeAdot0XXX * (
                mixDerivativeAdot00XX * (
                    mixDerivativeAdot000X * l0000Hash.x
                  + mixDerivativeBdot000X * l0001Hash.x
                )
                + mixDerivativeBdot00XX * (
                    mixDerivativeAdot001X * l0010Hash.x
                  + mixDerivativeBdot001X * l0011Hash.x
                )
            )
            + mixDerivativeBdot0XXX * (
                mixDerivativeAdot01XX * (
                    mixDerivativeAdot010X * l0100Hash.x
                  + mixDerivativeBdot010X * l0101Hash.x
                )
                + mixDerivativeBdot01XX * (
                    mixDerivativeAdot011X * l0110Hash.x
                  + mixDerivativeBdot011X * l0111Hash.x
                )
            )
        )
        + mixDerivativeBdotXXXX * (
            mixDerivativeAdot1XXX * (
                mixDerivativeAdot10XX * (
                    mixDerivativeAdot100X * l1000Hash.x
                  + mixDerivativeBdot100X * l1001Hash.x
                )
                + mixDerivativeBdot10XX * (
                    mixDerivativeAdot101X * l1010Hash.x
                  + mixDerivativeBdot101X * l1011Hash.x
                )
            )
            + mixDerivativeBdot1XXX * (
                mixDerivativeAdot11XX * (
                    mixDerivativeAdot110X * l1100Hash.x
                  + mixDerivativeBdot110X * l1101Hash.x
                )
                + mixDerivativeBdot11XX * (
                    mixDerivativeAdot111X * l1110Hash.x
                  + mixDerivativeBdot111X * l1111Hash.x
                )
            )
        )
        + mixDerivativeT(mixdot0XXX, mixdot1XXX, tMapped.x) * tMappedDerivative.x,

        // dz/dty
        mixDerivativeAdotXXXX * (
            mixDerivativeAdot0XXX * (
                mixDerivativeAdot00XX * (
                    mixDerivativeAdot000X * l0000Hash.y
                  + mixDerivativeBdot000X * l0001Hash.y
                )
                + mixDerivativeBdot00XX * (
                    mixDerivativeAdot001X * l0010Hash.y
                  + mixDerivativeBdot001X * l0011Hash.y
                )
            )
            + mixDerivativeBdot0XXX * (
                mixDerivativeAdot01XX * (
                    mixDerivativeAdot010X * l0100Hash.y
                  + mixDerivativeBdot010X * l0101Hash.y
                )
                + mixDerivativeBdot01XX * (
                    mixDerivativeAdot011X * l0110Hash.y
                  + mixDerivativeBdot011X * l0111Hash.y
                )
            )
            + mixDerivativeT(mixdot00XX, mixdot01XX, tMapped.y) * tMappedDerivative.y
        )
        + mixDerivativeBdotXXXX * (
            mixDerivativeAdot1XXX * (
                mixDerivativeAdot10XX * (
                    mixDerivativeAdot100X * l1000Hash.y
                  + mixDerivativeBdot100X * l1001Hash.y
                )
                + mixDerivativeBdot10XX * (
                    mixDerivativeAdot101X * l1010Hash.y
                  + mixDerivativeBdot101X * l1011Hash.y
                )
            )
            + mixDerivativeBdot1XXX * (
                mixDerivativeAdot11XX * (
                    mixDerivativeAdot110X * l1100Hash.y
                  + mixDerivativeBdot110X * l1101Hash.y
                )
                + mixDerivativeBdot11XX * (
                    mixDerivativeAdot111X * l1110Hash.y
                  + mixDerivativeBdot111X * l1111Hash.y
                )
            )
            + mixDerivativeT(mixdot10XX, mixdot11XX, tMapped.y) * tMappedDerivative.y
        ),

        // dz/dtz
        mixDerivativeAdotXXXX * (
            mixDerivativeAdot0XXX * (
                mixDerivativeAdot00XX * (
                    mixDerivativeAdot000X * l0000Hash.z
                  + mixDerivativeBdot000X * l0001Hash.z
                )
                + mixDerivativeBdot00XX * (
                    mixDerivativeAdot001X * l0010Hash.z
                  + mixDerivativeBdot001X * l0011Hash.z
                )
                + mixDerivativeT(mixdot000X, mixdot001X, tMapped.z) * tMappedDerivative.z
            )
            + mixDerivativeBdot0XXX * (
                mixDerivativeAdot01XX * (
                    mixDerivativeAdot010X * l0100Hash.z
                  + mixDerivativeBdot010X * l0101Hash.z
                )
                + mixDerivativeBdot01XX * (
                    mixDerivativeAdot011X * l0110Hash.z
                  + mixDerivativeBdot011X * l0111Hash.z
                )
                + mixDerivativeT(mixdot010X, mixdot011X, tMapped.z) * tMappedDerivative.z
            )
        )
        + mixDerivativeBdotXXXX * (
            mixDerivativeAdot1XXX * (
                mixDerivativeAdot10XX * (
                    mixDerivativeAdot100X * l1000Hash.z
                  + mixDerivativeBdot100X * l1001Hash.z
                )
                + mixDerivativeBdot10XX * (
                    mixDerivativeAdot101X * l1010Hash.z
                  + mixDerivativeBdot101X * l1011Hash.z
                )
                + mixDerivativeT(mixdot100X, mixdot101X, tMapped.z) * tMappedDerivative.z
            )
            + mixDerivativeBdot1XXX * (
                mixDerivativeAdot11XX * (
                    mixDerivativeAdot110X * l1100Hash.z
                  + mixDerivativeBdot110X * l1101Hash.z
                )
                + mixDerivativeBdot11XX * (
                    mixDerivativeAdot111X * l1110Hash.z
                  + mixDerivativeBdot111X * l1111Hash.z
                )
                + mixDerivativeT(mixdot110X, mixdot111X, tMapped.z) * tMappedDerivative.z
            )
        ),

        // dz/dtw
        mixDerivativeAdotXXXX * (
            mixDerivativeAdot0XXX * (
                mixDerivativeAdot00XX * (
                    mixDerivativeAdot000X * l0000Hash.w
                  + mixDerivativeBdot000X * l0001Hash.w
                  + mixDerivativeT(dot0000, dot0001, tMapped.w) * tMappedDerivative.w
                )
                + mixDerivativeBdot00XX * (
                    mixDerivativeAdot001X * l0010Hash.w
                  + mixDerivativeBdot001X * l0011Hash.w
                  + mixDerivativeT(dot0010, dot0011, tMapped.w) * tMappedDerivative.w
                )
            )
            + mixDerivativeBdot0XXX * (
                mixDerivativeAdot01XX * (
                    mixDerivativeAdot010X * l0100Hash.w
                  + mixDerivativeBdot010X * l0101Hash.w
                  + mixDerivativeT(dot0100, dot0101, tMapped.w) * tMappedDerivative.w
                )
                + mixDerivativeBdot01XX * (
                    mixDerivativeAdot011X * l0110Hash.w
                  + mixDerivativeBdot011X * l0111Hash.w
                  + mixDerivativeT(dot0110, dot0111, tMapped.w) * tMappedDerivative.w
                )
            )
        )
        + mixDerivativeBdotXXXX * (
            mixDerivativeAdot1XXX * (
                mixDerivativeAdot10XX * (
                    mixDerivativeAdot100X * l1000Hash.w
                  + mixDerivativeBdot100X * l1001Hash.w
                  + mixDerivativeT(dot1000, dot1001, tMapped.w) * tMappedDerivative.w
                )
                + mixDerivativeBdot10XX * (
                    mixDerivativeAdot101X * l1010Hash.w
                  + mixDerivativeBdot101X * l1011Hash.w
                  + mixDerivativeT(dot1010, dot1011, tMapped.w) * tMappedDerivative.w
                )
            )
            + mixDerivativeBdot1XXX * (
                mixDerivativeAdot11XX * (
                    mixDerivativeAdot110X * l1100Hash.w
                  + mixDerivativeBdot110X * l1101Hash.w
                  + mixDerivativeT(dot1100, dot1101, tMapped.w) * tMappedDerivative.w
                )
                + mixDerivativeBdot11XX * (
                    mixDerivativeAdot111X * l1110Hash.w
                  + mixDerivativeBdot111X * l1111Hash.w
                  + mixDerivativeT(dot1110, dot1111, tMapped.w) * tMappedDerivative.w
                )
            )
        )
    );
}

float perlinNoiseUnsigned1(float x) { return (perlinNoiseSigned1(x) + 1.0) * 0.5; }
float perlinNoiseUnsigned2(vec2 v)  { return (perlinNoiseSigned2(v) + 1.0) * 0.5; }
float perlinNoiseUnsigned3(vec3 v)  { return (perlinNoiseSigned3(v) + 1.0) * 0.5; }
float perlinNoiseUnsigned4(vec4 v)  { return (perlinNoiseSigned4(v) + 1.0) * 0.5; }
float perlinNoiseUnsignedGradient1(float x) { return perlinNoiseSignedGradient1(x) * 0.5; }
vec2  perlinNoiseUnsignedGradient2(vec2 v)  { return perlinNoiseSignedGradient2(v) * 0.5; }
vec3  perlinNoiseUnsignedGradient3(vec3 v)  { return perlinNoiseSignedGradient3(v) * 0.5; }
vec4  perlinNoiseUnsignedGradient4(vec4 v)  { return perlinNoiseSignedGradient4(v) * 0.5; }

/*
 * Macros to compute fractal versions of noise functions with the following parameters:
 *
 * tyAssignTo: The type representing the dimension of the noise function (typically `float`, but also `float, `vec2`, `vec3`, `vec4` for noise function derivatives);
 * tyInput: The type of `exprV` and consequently the type of the single argument of `identNoiseFunction`;
 * identAssignTo: The name of the variable of type `tyAssignTo` to assign the result to;
 * identNoiseFunction: The name of the base n-dimensional noise function;
 * exprV: The input vector;
 * exprLacunarity: The modifier to multiply `exprV` each iteration;
 * exprGain: The modifier to multiply `exprV` each iteration;
 * exprLayers: The number of layers of the noise to sum up;
 */
#define FRACTALIFY(tyAssignTo, tyInput, identAssignTo, identNoiseFunction, exprV, exprLacunarity, exprGain, exprLayers) ;\
tyAssignTo identAssignTo;                                                          \
do {                                                                               \
    tyAssignTo result = tyAssignTo(0.0);                                           \
    tyInput    inMultiplier = tyInput(1.0);                                        \
    tyAssignTo outMultiplier = tyAssignTo(1.0);                                    \
    tyAssignTo range = tyAssignTo(0.0);                                            \
                                                                                   \
    for (uint i = 0; i < uint(exprLayers); i++) {                                  \
        result += identNoiseFunction((exprV) * inMultiplier) * outMultiplier;      \
        range += outMultiplier;                                                    \
        inMultiplier *= tyInput(exprLacunarity);                                   \
        outMultiplier *= tyAssignTo(exprGain);                                     \
    }                                                                              \
                                                                                   \
    tyAssignTo normalizedResult = result / range;                                  \
    identAssignTo = normalizedResult;                                              \
} while(false);

// A version of `FRACTALIFY`, where `exprGain = 1 / exprLacunarity`
#define FRACTALIFY_PINK(tyAssignTo, tyInput, identAssignTo, identNoiseFunction, exprV, exprLacunarity, exprLayers) \
    FRACTALIFY(tyAssignTo, tyInput, identAssignTo, identNoiseFunction, exprV, exprLacunarity, (1.0 / exprLacunarity), exprLayers)

// A version of `FRACTALIFY`, where `exprLacunarity = 2` and `exprGain = 0.5`
#define FRACTALIFY_BROWN(tyAssignTo, tyInput, identAssignTo, identNoiseFunction, exprV, exprLayers) \
    FRACTALIFY_PINK(tyAssignTo, tyInput, identAssignTo, identNoiseFunction, exprV, 2.0, exprLayers)
