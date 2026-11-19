// Safe WGSL shaders for ConnectomeAtlas v0.5.0:
// 1. Anisotropic Banks/Hege illuminated braided white-matter streamlines with calm tonic laminar wave
// 2. Fresnel rim-lit continuous cortical surface heatmap with CSG surgical octant/hemisphere cutaway
struct ViewUniforms {
    view_proj: mat4x4<f32>,
    camera_eye: vec4<f32>,
    tonic_phase: f32,
    shell_opacity: f32,
    cutaway_offset: f32,
    cutaway_mode: u32,
};

@group(0) @binding(0)
var<uniform> uniforms: ViewUniforms;

struct StreamlineVertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) tangent: vec3<f32>,
    @location(3) color: vec4<f32>,
    @location(4) arclength_u: f32,
};

struct StreamlineVertexOutput {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
    @location(1) world_normal: vec3<f32>,
    @location(2) world_tangent: vec3<f32>,
    @location(3) color: vec4<f32>,
    @location(4) arclength_u: f32,
};

@vertex
fn vs_main(in: StreamlineVertexInput) -> StreamlineVertexOutput {
    var out: StreamlineVertexOutput;
    out.clip_pos = uniforms.view_proj * vec4<f32>(in.position, 1.0);
    out.world_pos = in.position;
    out.world_normal = normalize(in.normal);
    out.world_tangent = normalize(in.tangent);
    out.color = in.color;
    out.arclength_u = in.arclength_u;
    return out;
}

@fragment
fn fs_main(in: StreamlineVertexOutput) -> @location(0) vec4<f32> {
    let t = normalize(in.world_tangent);
    let n = normalize(in.world_normal);
    let l = normalize(vec3<f32>(0.48, 0.56, 0.67));
    let v = normalize(uniforms.camera_eye.xyz - in.world_pos);
    let lt = clamp(dot(l, t), -1.0, 1.0);
    let vt = clamp(dot(v, t), -1.0, 1.0);
    let sin_lt = sqrt(max(0.0, 1.0 - lt * lt));
    let sin_vt = sqrt(max(0.0, 1.0 - vt * vt));
    let aniso_diff = 0.36 + 0.64 * sin_lt;
    let aniso_spec = pow(max(0.0, sin_lt * sin_vt - lt * vt), 24.0);
    let cyl_shade = 0.65 + 0.35 * max(0.0, dot(n, l));
    let cam_dist = length(uniforms.camera_eye.xyz - in.world_pos);
    let depth_cue = clamp(1.14 - (cam_dist - 125.0) / 360.0, 0.64, 1.12);
    let laminar_wave = 0.87 + 0.13 * sin((in.arclength_u - uniforms.tonic_phase) * 6.2831853);
    let lit = in.color.rgb * (aniso_diff * cyl_shade * laminar_wave) * depth_cue
            + vec3<f32>(0.96, 0.98, 1.0) * (aniso_spec * 0.34 * depth_cue);
    return vec4<f32>(lit, in.color.a);
}

@fragment
fn fs_fresnel_cortex(in: StreamlineVertexOutput) -> @location(0) vec4<f32> {
    var incision_glow: f32 = 0.0;
    if (uniforms.cutaway_mode == 1u) {
        let cx = 2.0 + uniforms.cutaway_offset * 0.45;
        let cy = -6.0 - uniforms.cutaway_offset * 0.55;
        let cz = 0.0 - uniforms.cutaway_offset * 0.45;
        if (in.world_pos.x < cx && in.world_pos.y > cy && in.world_pos.z > cz) {
            discard;
        }
        let d_corner = max(max(cx - in.world_pos.x, in.world_pos.y - cy), in.world_pos.z - cz);
        if (abs(d_corner) < 2.2) {
            incision_glow = 1.0 - abs(d_corner) / 2.2;
        }
    } else if (uniforms.cutaway_mode == 2u) {
        if (in.world_pos.x < uniforms.cutaway_offset) {
            discard;
        }
    } else if (uniforms.cutaway_mode == 3u) {
        if (in.world_pos.z > (16.0 - uniforms.cutaway_offset * 0.65)) {
            discard;
        }
    }
    let n = normalize(in.world_normal);
    let v = normalize(uniforms.camera_eye.xyz - in.world_pos);
    let fresnel = pow(1.0 - clamp(abs(dot(n, v)), 0.0, 1.0), 2.10);
    let col = mix(in.color.rgb + vec3<f32>(0.85, 0.95, 1.0) * (fresnel * 0.36), vec3<f32>(0.35, 0.92, 1.0), incision_glow * 0.8);
    let alpha = clamp(uniforms.shell_opacity * (0.34 + 1.72 * fresnel) + incision_glow * 0.52, 0.02, 0.92);
    return vec4<f32>(col, alpha);
}
