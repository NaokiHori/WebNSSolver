#version 300 es

precision mediump float;

const float PI = 3.14159265358979;

uniform sampler2D u_texture;
uniform int u_scalar_type;

in vec2 v_texture_coordinates;

out vec4 frag_color;

void main(void) {
  float value = texture(u_texture, v_texture_coordinates).r;
  if (0 == u_scalar_type) {
    float r = value < 0.5 ? 1. : 2. - 2. * value;
    float g = value < 0.5 ? 2. * value : 2. - 2. * value;
    float b = value < 0.5 ? 2. * value : 1.;
    frag_color = vec4(r, g, b, 1.);
  } else if (1 == u_scalar_type) {
    float r = value < 0.5 ? 1. - 2. * value : 0.;
    float g = value < 0.5 ? 1. - 2. * value : 2. * value - 1.;
    float b = value < 0.5 ? 0. : 2. * value - 1.;
    frag_color = vec4(r, g, b, 1.);
  } else if (2 == u_scalar_type) {
    float g = value < 0.5 ? 1. - 2. * value : 0.;
    float b = value < 0.5 ? 1. - 2. * value : 2. * value - 1.;
    float r = value < 0.5 ? 0. : 2. * value - 1.;
    frag_color = vec4(r, g, b, 1.);
  } else {
    frag_color = vec4(vec3(value), 1.);
  }
}

