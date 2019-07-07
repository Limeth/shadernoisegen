#extension GL_ARB_derivative_control : enable
#include "noise_configuration.glsl"
SHADERNOISEGEN(quintic)
#include "noise.glsl"

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Normalized pixel coordinates (from 0 to 1)
  vec2 uv = fragCoord / iResolution.xy;
  float u = uv.x;
  float v = uv.y;
  vec4 color = vec4(vec3(0.0), 1.0);

  /* float value; FRACTALIFY_BROWN(value, valueNoiseUnsigned4, vec4(uv * 10, iTime * 0.2, 0), 10) */
  /* float value; FRACTALIFY_BROWN(value, perlinNoiseUnsigned4, vec4(uv * 10, iTime * 1.0, (uv.x + uv.y) * 0.0), 10) */
  uint layers = 7;
  vec4 arg = vec4(0, uv * 10, iTime * 0.2);
  /* float value = valueNoiseUnsigned3(arg); */
  FRACTALIFY_BROWN(float, value, valueNoiseUnsigned4, arg, layers);
  /* vec3 deriv = valueNoiseUnsignedGradient3(arg); */
  FRACTALIFY_BROWN(vec4, deriv, valueNoiseUnsignedGradient4, arg, layers);

  /* value = 0; */
  /* deriv = vec2(0); */

  deriv.xy = mod((fragCoord.x + fragCoord.y) / 2, 2) == 0 ? abs(deriv.xy) : deriv.xy;
  deriv.z = mod((fragCoord.x + fragCoord.y) / 4, 4) == 0 ? deriv.z : 0;
  deriv.w = mod((fragCoord.x - fragCoord.y) / 4, 4) == 0 ? deriv.w : 0;
  /* deriv.x = 0; */
  /* deriv = vec2(dFdxFine(value), dFdyFine(value)) * 50; */

  /* deriv = vec2(0); */
  /* float deriv = valueNoiseUnsignedDerivative1(arg); */
  /* float deriv = valueNoiseUnsignedDerivativeAlt1(arg); */
  color.r = value;
  /* color.gb += deriv.xy; */
  /* color.rgb += vec3(deriv.z) * 2.0; */
  color.rgb += vec3(deriv.w) * 2.0;

  // Output to screen
  fragColor = color;
}

out vec4 fragColor;
void main() { mainImage(fragColor, gl_FragCoord.xy); }
