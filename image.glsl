#extension GL_ARB_derivative_control : enable
#include "noise_configuration.glsl"
SHADERNOISEGEN(linear)
#include "noise.glsl"

float turbulentValueNoise3(vec3 x) {
    return abs(valueNoiseSigned3(x));
}

float turbulentPerlinNoise3(vec3 x) {
    return abs(perlinNoiseSigned3(x));
}

vec3 turbulentPerlinNoiseGradient3(vec3 x) {
    float s = sign(perlinNoiseSigned3(x));
    return perlinNoiseSignedGradient3(x) * s;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Normalized pixel coordinates (from 0 to 1)
  vec2 uv = fragCoord / max(iResolution.x, iResolution.y);
  float u = uv.x;
  float v = uv.y;
  vec4 color = vec4(vec3(0.0), 1.0);

  /* float value; FRACTALIFY_BROWN(value, valueNoiseUnsigned4, vec4(uv * 10, iTime * 0.2, 0), 10) */
  /* float value; FRACTALIFY_BROWN(value, perlinNoiseUnsigned4, vec4(uv * 10, iTime * 1.0, (uv.x + uv.y) * 0.0), 10) */
  uint layers = 16;
  /* float arg = u * 20 + 20; */
  /* vec2 arg = vec2(0, u * 10); */
  /* vec2 arg = uv * 10 + vec2(0, 0); */
  /* vec4 arg = vec4(uv * 10, iTime * 0.2, 0).zwxy; */
  vec3 arg = vec3(uv * 1, iTime * 0.2);
  /* vec3 arg = vec3(u + iTime * 0.2, 0, 0); */
  const float GOLDEN_RATIO = 1.6180339887498948482;
  /* float lacunarity = float(GOLDEN_RATIO); */
  float lacunarity = GOLDEN_RATIO;
  float gain = 1.0 / lacunarity;
  /* float value = valueNoiseUnsigned3(arg); */
  /* FRACTALIFY(float, float, value, perlinNoiseUnsigned1, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, float, grad, perlinNoiseUnsignedGradient1, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, float, value, valueNoiseUnsigned1, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, float, grad, valueNoiseUnsignedGradient1, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, vec3, value, turbulentPerlinNoise3, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(vec3, vec3, grad, turbulentPerlinNoiseGradient3, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, vec2, value, worleyNoiseUnsignedSingle2, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, vec3, value, worleyNoiseUnsignedDouble3, arg, lacunarity, gain, layers); */
  FRACTALIFY(vec3,  vec3, grad, perlinNoiseUnsignedGradient3, arg, lacunarity, gain, layers);
  /* FRACTALIFY(float, vec2, value, worleyNoiseUnsignedSingle2, arg, lacunarity, gain, layers); */
  /* value = pow(atan(value), 2); */
  /* FRACTALIFY(vec3, vec3, grad, worleyNoiseUnsignedGradient1, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(float, vec3, value, valueNoiseUnsigned3, arg, lacunarity, gain, layers); */
  /* FRACTALIFY(vec3, vec3, grad, valueNoiseUnsignedGradient3, arg, lacunarity, gain, layers); */

  /* vec2 lowResCoord = floor(fragCoord / 4); */
  /* /1* grad.xyzw = grad.zwxy; *1/ */
  /* grad.x = mod(lowResCoord.x + lowResCoord.y, 4) == 0 ? abs(grad.x) : grad.x; */
  /* grad.y = mod(lowResCoord.x + lowResCoord.y + 2, 4) == 0 ? abs(grad.y) : grad.y; */
  float y = fragCoord.y / iResolution.y;
  float lineThickness = 0.01;
  /* grad *= 0.2;// * (sin(iTime) / 2 + 0.5); */
  /* fragColor.rgb = vec3(value); */
  /* if (uv.x > 0.5) { */
  /*     fragColor.rg += grad.xy; */
  /* } */

  /* fragColor.rgb = vec3(1 - step(value, y)); */
  /* fragColor.rgb = vec3(value); */
  /* fragColor.r = (abs(grad - y + 0.5) < lineThickness / 2.0) ? 1 : 0; */
  /* fragColor.rgb += vec3(abs(y - 0.5) < lineThickness / 12.0) / 2; */

  /* vec3 grad = valueNoiseUnsignedGradient3(arg); */
  /* FRACTALIFY(vec3, vec3, grad, valueNoiseUnsignedGradient3, arg, lacunarity, gain, layers); */

  /* value = 0; */
  /* grad = vec2(0); */

  /* grad.xy = mod((fragCoord.x + fragCoord.y) / 2, 2) == 0 ? abs(grad.xy) : grad.xy; */
  /* grad.z = mod((fragCoord.x + fragCoord.y) / 4, 4) == 0 ? grad.z : 0; */
  /* grad.w = mod((fragCoord.x - fragCoord.y) / 4, 4) == 0 ? grad.w : 0; */
  /* grad.x = 0; */

  /* grad = vec2(0); */
  /* float grad = valueNoiseUnsignedgradative1(arg); */
  /* float grad = valueNoiseUnsignedgradativeAlt1(arg); */
  /* fragColor.rgb = vec3(value); */
  /* fragColor.rgb = vec3(grad.z); */
  /* fragColor.gb += grad.xy; */
  /* if (fragCoord.y / iResolution.y < 0.5) { */
  /*     fragColor.g += abs(grad.x); */
  /* } else { */
  /*     fragColor.b += abs(grad.y); */
  /* } */
  /* fragColor.rgb += vec3(grad.z) * 2.0; */
  /* fragColor.rgb += vec3(grad.w) * 2.0; */

  vec2 texCoord = fragCoord / iResolution.xy;
  texCoord += -grad.xy * 0.150;

  /* fragColor.rgb = grad; */


  // Output to screen
  /* if ((fragCoord.x / iResolution.x) > 0.5) { */
      fragColor = texture(background, texCoord);
  /* } else { */
      /* fragColor += vec4(vec3(value), 0); */
  /* } */

  /* fragColor.r += 0.5 - step(grad + 0.5, y) / 2.0; */
}

out vec4 fragColor;
void main() { mainImage(fragColor, gl_FragCoord.xy); }
