#include "noise.glsl"

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Normalized pixel coordinates (from 0 to 1)
  vec2 uv = fragCoord / iResolution.xy;
  // Time varying pixel color
  vec4 color = vec4(vec3(0.0), 1.0);
  /* color.rgb = vec3(valueNoiseUnsigned4(vec4(uv * 100.0, iTime, iTime * (uv.x - uv.y)))); */
  /* color.rgb = vec3(perlinNoiseUnsigned3(vec3(uv * 10, iTime))); */
  float lacunarity = uv.x;
  lacunarity = floor(lacunarity * 10) / 10;
  lacunarity *= 3;
  lacunarity = pow(2, lacunarity);
  float gain = uv.y;
  gain = floor(gain * 10) / 10;
  gain *= 3;
  gain = pow(2, gain);
  gain = 1 / gain;
  float noise;
  /* FRACTALIFY(noise, valueNoiseUnsigned3, vec3(uv * 10, iTime * 1.0), lacunarity, gain, 10) */
  /* FRACTALIFY_PINK(noise, valueNoiseUnsigned3, vec3(uv * 10, iTime * 1.0), lacunarity, 10) */
  FRACTALIFY_BROWN(noise, valueNoiseUnsigned3, vec3(uv * 10, iTime * 0.2), 10)

  color.rgb = vec3(noise);
  // Output to screen
  fragColor = color;
}

out vec4 fragColor;
void main() { mainImage(fragColor, gl_FragCoord.xy); }
