const float TAU = 6.28318530718;
const float PI  = 3.14159265359;

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

// TODO: Figure out how to generate a deterministically random 4D unit vector
/* vec4 ihash3AsUnitVec3(vec4 v); */

float smoothStep1(float t) { return t * t * (3 - 2 * t); }
vec2  smoothStep2(vec2 t)  { return t * t * (3 - 2 * t); }
vec3  smoothStep3(vec3 t)  { return t * t * (3 - 2 * t); }
vec4  smoothStep4(vec4 t)  { return t * t * (3 - 2 * t); }

float quintic1(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec2  quintic2(vec2  t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec3  quintic3(vec3  t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec4  quintic4(vec4  t) { return t * t * t * (t * (t * 6 - 15) + 10); }

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
    float tMapped = smoothStep1(t);

    return mix(l0Hash, l1Hash, tMapped);
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
    vec2 tMapped = smoothStep2(t);

    return mix(
        mix(l00Hash, l01Hash, tMapped.y),
        mix(l10Hash, l11Hash, tMapped.y),
        tMapped.x
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
    vec3 tMapped = smoothStep3(t);

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
    vec4 tMapped = smoothStep4(t);

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

float valueNoiseSigned1(float x) { return valueNoiseUnsigned1(x) * 2.0 - 1.0; }
float valueNoiseSigned2(vec2 v) { return valueNoiseUnsigned2(v) * 2.0 - 1.0; }
float valueNoiseSigned3(vec3 v) { return valueNoiseUnsigned3(v) * 2.0 - 1.0; }
float valueNoiseSigned4(vec4 v) { return valueNoiseUnsigned4(v) * 2.0 - 1.0; }

float perlinNoiseSigned1(float x) {
    // Lattice point integer coordinates
    int l0 = int(floor(x));
    int l1 = l0 + 1;
    // Deterministically random values for lattice points
    float l0Hash = ihash1AsUnitVec1(l0);
    float l1Hash = ihash1AsUnitVec1(l1);
    // Vectors to X from each lattice point
    float v0 = fract(x);
    float v1 = float(v0.x - 1);
    // Interpolation
    float t = v0;
    float tMapped = smoothStep1(t);

    return mix(dot(l0Hash, v0), dot(l1Hash, v1), tMapped);
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
    vec2 tMapped = quintic2(t);

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
    vec3 tMapped = quintic3(t);

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

float perlinNoiseUnsigned1(float x) { return (perlinNoiseSigned1(x) + 1.0) * 0.5; }
float perlinNoiseUnsigned2(vec2 v)  { return (perlinNoiseSigned2(v) + 1.0) * 0.5; }
float perlinNoiseUnsigned3(vec3 v)  { return (perlinNoiseSigned3(v) + 1.0) * 0.5; }

/*
 * Implementations of fractal sums with the following parameters:
 * x or v:     Input vector;
 * lacunarity: Input multiplier for each layer, typically > 1;
 * gain:       Output multiplier for each layer, typically in (0; 1);
 */

float valueNoiseFractalUnsigned1(float x, float lacunarity, float gain, uint layers) {
    float result = 0.0;
    float inMultiplier = 1.0;
    float outMultiplier = 1.0;
    float range = 0.0;

    for (uint i = 0; i < layers; i++) {
        result += valueNoiseUnsigned1(x * inMultiplier) * outMultiplier;
        range += outMultiplier;
        inMultiplier *= lacunarity;
        outMultiplier *= gain;
    }

    float normalizedResult = result / range;

    return normalizedResult;
}

/*
 * Macros to compute fractal versions of noise functions with the following parameters:
 *
 * identAssignTo: The name of the variable of type `float` to assign the result to;
 * identNoiseFunction: The name of the base n-dimensional noise function;
 * exprV: The input vector;
 * exprLacunarity: The modifier to multiply `exprV` each iteration;
 * exprGain: The modifier to multiply `exprV` each iteration;
 * exprLayers: The number of layers of the noise to sum up;
 */
#define FRACTALIFY(identAssignTo, identNoiseFunction, exprV, exprLacunarity, exprGain, exprLayers) do { \
    float result = 0.0;                                                            \
    float inMultiplier = 1.0;                                                      \
    float outMultiplier = 1.0;                                                     \
    float range = 0.0;                                                             \
                                                                                   \
    for (uint i = 0; i < uint(exprLayers); i++) {                                  \
        result += identNoiseFunction((exprV) * inMultiplier) * outMultiplier;      \
        range += outMultiplier;                                                    \
        inMultiplier *= float(exprLacunarity);                                     \
        outMultiplier *= float(exprGain);                                          \
    }                                                                              \
                                                                                   \
    float normalizedResult = result / range;                                       \
    identAssignTo = normalizedResult;                                              \
} while(false);

// A version of `FRACTALIFY`, where `exprGain = 1 / exprLacunarity`
#define FRACTALIFY_PINK(identAssignTo, identNoiseFunction, exprV, exprLacunarity, exprLayers) \
    FRACTALIFY(identAssignTo, identNoiseFunction, exprV, exprLacunarity, (1.0 / exprLacunarity), exprLayers)

// A version of `FRACTALIFY`, where `exprLacunarity = 2` and `exprGain = 0.5`
#define FRACTALIFY_BROWN(identAssignTo, identNoiseFunction, exprV, exprLayers) \
    FRACTALIFY_PINK(identAssignTo, identNoiseFunction, exprV, 2.0, exprLayers)
