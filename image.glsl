#extension GL_ARB_derivative_control : enable
#include "noise_configuration.glsl"
SHADERNOISEGEN(quintic)
#include "noise.glsl"

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Normalized pixel coordinates (from 0 to 1)
  vec2 uv = fragCoord / max(iResolution.x, iResolution.y);
  float u = uv.x;
  float v = uv.y;
  vec4 color = vec4(vec3(0.0), 1.0);

  /* float value; FRACTALIFY_BROWN(value, valueNoiseUnsigned4, vec4(uv * 10, iTime * 0.2, 0), 10) */
  /* float value; FRACTALIFY_BROWN(value, perlinNoiseUnsigned4, vec4(uv * 10, iTime * 1.0, (uv.x + uv.y) * 0.0), 10) */
  uint layers = 16;
  vec3 arg = vec3(uv * 10, iTime * 0.5);
  /* vec3 arg = vec3(u + iTime * 0.2, 0, 0); */
  const float GOLDEN_RATIO = 1.6180339887498948482;
  vec3 lacunarity = vec3(GOLDEN_RATIO);

  vec3 gain = 1.0 / lacunarity;
  /* float value = valueNoiseUnsigned3(arg); */
  float value;
  vec3 grad;

  if (fragCoord.x / iResolution.x > 1.0) {
      FRACTALIFY(float, vec3, valuePerlin, perlinNoiseUnsigned3, arg, lacunarity, gain, layers);
      /* FRACTALIFY(vec3, vec3, gradPerlin, perlinNoiseUnsignedGradient3, arg, lacunarity, gain, layers); */
      value = valuePerlin;
      grad = vec3(0);

  } else {
      FRACTALIFY(float, vec3, valueValue,  valueNoiseUnsigned3, arg, lacunarity, gain, layers);
      FRACTALIFY(vec3, vec3, gradValue, valueNoiseUnsignedGradient3, arg, lacunarity, gain, layers);
      value = valueValue;
      grad = gradValue;
  }

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
  /* color.rgb = vec3(value); */
  /* color.rgb = vec3(grad.z); */
  /* color.gb += grad.xy; */
  /* if (fragCoord.y / iResolution.y < 0.5) { */
  /*     color.g += abs(grad.x); */
  /* } else { */
  /*     color.b += abs(grad.y); */
  /* } */
  /* color.rgb += vec3(grad.z) * 2.0; */
  /* color.rgb += vec3(grad.w) * 2.0; */

  vec2 texCoord = fragCoord / iResolution.xy;
  texCoord += grad.xy * 0.025;

  /* color.rgb = grad; */


  // Output to screen
  /* if (mod((fragCoord.x + fragCoord.y) / 2, 2) != 0) { */
      fragColor = texture(background, texCoord);
  /* } else { */
      /* fragColor += color * 1; */
  /* } */

  /* fragColor.rgb = vec3(step(value, (fragCoord / iResolution.xy).y)); */
}

out vec4 fragColor;
void main() { mainImage(fragColor, gl_FragCoord.xy); }
