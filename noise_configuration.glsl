const float TAU = 6.28318530718;
const float PI  = 3.14159265359;

#define SHADERNOISEGEN(interpolationFunction)                                            \
float noiseInterpolation(float t) { return interpolationFunction(t); }                   \
vec2  noiseInterpolation(vec2  t) { return interpolationFunction(t); }                   \
vec3  noiseInterpolation(vec3  t) { return interpolationFunction(t); }                   \
vec4  noiseInterpolation(vec4  t) { return interpolationFunction(t); }                   \
float noiseInterpolationGradient(float t) { return interpolationFunction##Gradient(t); } \
vec2  noiseInterpolationGradient(vec2  t) { return interpolationFunction##Gradient(t); } \
vec3  noiseInterpolationGradient(vec3  t) { return interpolationFunction##Gradient(t); } \
vec4  noiseInterpolationGradient(vec4  t) { return interpolationFunction##Gradient(t); } \

float smoothStep(float t) { return t * t * (3 - 2 * t); }
vec2  smoothStep(vec2 t)  { return t * t * (3 - 2 * t); }
vec3  smoothStep(vec3 t)  { return t * t * (3 - 2 * t); }
vec4  smoothStep(vec4 t)  { return t * t * (3 - 2 * t); }

/*
 * Derivation:
 * smoothStep(t) = t * t * (3 - 2 * t)
 * smoothStep(t) = 3t^2 - 2t^3
 * smoothStep'(t) = 6t - 6t^2
 * smoothStep'(t) = 6t * (1 - t)
 */
float smoothStepGradient(float t) { return 6 * t * (1 - t); }
vec2  smoothStepGradient(vec2 t)  { return 6 * t * (1 - t); }
vec3  smoothStepGradient(vec3 t)  { return 6 * t * (1 - t); }
vec4  smoothStepGradient(vec4 t)  { return 6 * t * (1 - t); }

float quintic(float t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec2  quintic(vec2  t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec3  quintic(vec3  t) { return t * t * t * (t * (t * 6 - 15) + 10); }
vec4  quintic(vec4  t) { return t * t * t * (t * (t * 6 - 15) + 10); }

/*
 * Derivation:
 * quintic(t) = t * t * t * (t * (t * 6 - 15) + 10)
 * quintic(t) = 6t^5 - 15t^4 + 10t^3
 * quintic'(t) = 30t^4 - 60t^3 + 30t^2
 * quintic'(t) = 30t^2 * (t^2 - 2t + 1)
 * quintic'(t) = 30t^2 * (t - 1)^2
 */
float quinticGradient(float t) {
    float tMinusOne = t - 1;
    return 30 * t * t * tMinusOne * tMinusOne;
}

vec2 quinticGradient(vec2 t) {
    vec2 tMinusOne = t - 1;
    return 30 * t * t * tMinusOne * tMinusOne;
}

vec3 quinticGradient(vec3 t) {
    vec3 tMinusOne = t - 1;
    return 30 * t * t * tMinusOne * tMinusOne;
}

vec4 quinticGradient(vec4 t) {
    vec4 tMinusOne = t - 1;
    return 30 * t * t * tMinusOne * tMinusOne;
}
