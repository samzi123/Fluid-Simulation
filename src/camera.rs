use bevy::prelude::*;
use bevy::window::{PrimaryWindow};
use bevy::input::mouse::{MouseWheel, MouseMotion};
use crate::MyWorldCoords;

/// Used to help identify our main camera
#[derive(Component)]
pub struct MainCamera;

/// Tags an entity as capable of panning and orbiting.
#[derive(Component)]
pub struct PanOrbitCamera {
    /// The "focus point" to orbit around. It is automatically updated when panning the camera
    pub focus: Vec3,
    pub radius: f32,
    pub upside_down: bool,
}

impl Default for PanOrbitCamera {
    fn default() -> Self {
        PanOrbitCamera {
            focus: Vec3::ZERO,
            radius: 5.0,
            upside_down: false,
        }
    }
}

pub fn spawn_camera_and_lights(commands: &mut Commands) {
    let camera_transform = Transform::from_xyz(5.0, 7.5, 30.0).looking_at(Vec3::ZERO, Vec3::Y);
    let rotation_radius = camera_transform.translation.length();

    // spawn camera
    commands.spawn((
        Camera3dBundle {
        transform: camera_transform,
        ..default()
        }, 
        MainCamera, 
        PanOrbitCamera {
            radius: rotation_radius,
            ..Default::default()
        },
    ));

    // spawn lights
    commands.spawn((
        PointLightBundle {
            point_light: PointLight {
                // shadows_enabled: true,
                intensity: 10000.0,
                ..default()
            },
            transform: Transform::from_xyz(4.0, 8.0, 4.0).looking_at(Vec3::ZERO, Vec3::Y),
            ..default()
        },
    ));
    commands.spawn((
        PointLightBundle {
            point_light: PointLight {
                // shadows_enabled: true,
                intensity: 10000.0,
                ..default()
            },
            transform: Transform::from_xyz(-4.0, -8.0, -10.0).looking_at(Vec3::ZERO, Vec3::Y),
            ..default()
        },
    ));
}

/// Pan the camera with right mouse click, zoom with scroll wheel, orbit with left mouse click.
pub fn pan_orbit_camera(
    q_window: Query<&Window, With<PrimaryWindow>>,
    mut ev_motion: EventReader<MouseMotion>,
    mut ev_scroll: EventReader<MouseWheel>,
    input_mouse: Res<Input<MouseButton>>,
    mut query: Query<(&mut PanOrbitCamera, &mut Transform, &Projection)>,
) {
    // change input mapping for orbit and panning here
    let orbit_button = MouseButton::Left;
    let pan_button = MouseButton::Right;

    let mut pan = Vec2::ZERO;
    let mut rotation_move = Vec2::ZERO;
    let mut scroll = 0.0;
    let mut orbit_button_changed = false;

    if input_mouse.pressed(orbit_button) {
        for ev in ev_motion.read() {
            rotation_move += ev.delta * crate::MOUSE_ROTATE_SPEED;
        }
    } else if input_mouse.pressed(pan_button) {
        // Pan only if we're not rotating at the moment
        for ev in ev_motion.read() {
            pan += ev.delta * crate::MOUSE_PAN_SPEED;
        }
    }
    for ev in ev_scroll.read() {
        scroll += ev.y * crate::MOUSE_SCROLL_SPEED;
    }
    if input_mouse.just_released(orbit_button) || input_mouse.just_pressed(orbit_button) {
        orbit_button_changed = true;
    }

    for (mut pan_orbit, mut transform, projection) in query.iter_mut() {
        if orbit_button_changed {
            // only check for upside down when orbiting started or ended this frame
            // if the camera is "upside" down, panning horizontally would be inverted, so invert the input to make it correct
            let up = transform.rotation * Vec3::Y;
            pan_orbit.upside_down = up.y <= 0.0;
        }

        let mut any = false;
        if rotation_move.length_squared() > 0.0 {
            any = true;
            let window = get_primary_window_size(&q_window);
            let delta_x = {
                let delta = rotation_move.x / window.x * std::f32::consts::PI * 2.0;
                if pan_orbit.upside_down { -delta } else { delta }
            };
            let delta_y = rotation_move.y / window.y * std::f32::consts::PI;
            let yaw = Quat::from_rotation_y(-delta_x);
            let pitch = Quat::from_rotation_x(-delta_y);
            transform.rotation = yaw * transform.rotation; // rotate around global y axis
            transform.rotation = transform.rotation * pitch; // rotate around local x axis
        } else if pan.length_squared() > 0.0 {
            any = true;
            // make panning distance independent of resolution and FOV,
            let window = get_primary_window_size(&q_window);
            if let Projection::Perspective(projection) = projection {
                pan *= Vec2::new(projection.fov * projection.aspect_ratio, projection.fov) / window;
            }
            // translate by local axes
            let right = transform.rotation * Vec3::X * -pan.x;
            let up = transform.rotation * Vec3::Y * pan.y;
            // make panning proportional to distance away from focus point
            let translation = (right + up) * pan_orbit.radius;
            pan_orbit.focus += translation;
        } else if scroll.abs() > 0.0 {
            any = true;
            pan_orbit.radius -= scroll * pan_orbit.radius * 0.2;
            // dont allow zoom to reach zero or you get stuck
            pan_orbit.radius = f32::max(pan_orbit.radius, 0.05);
        }

        if any {
            // emulating parent/child to make the yaw/y-axis rotation behave like a turntable
            // parent = x and y rotation
            // child = z-offset
            let rot_matrix = Mat3::from_quat(transform.rotation);
            transform.translation = pan_orbit.focus + rot_matrix.mul_vec3(Vec3::new(0.0, 0.0, pan_orbit.radius));
        }
    }

    // consume any remaining events, so they don't pile up if we don't need them
    // (and also to avoid Bevy warning us about not checking events every frame update)
    ev_motion.clear();
}

// Returns the mouse position in world coordinates
pub fn get_mouse_world_position(
    mycoords: &mut ResMut<MyWorldCoords>,
    // query to get the window (so we can read the current cursor position)
    window_state: &Window,
    // query to get camera transform
    q_camera: &Query<(&Camera, &GlobalTransform), With<MainCamera>>,
) -> Option<Vec3> {
    // get the camera info and transform
    // assuming there is exactly one main camera entity, so Query::single() is OK
    let (camera, camera_transform) = q_camera.single();

    // check if the cursor is inside the window and get its position
    // then, ask bevy to convert into world coordinates, and truncate to discard Z
    if let Some(world_position) = window_state.cursor_position()
        .and_then(|cursor| camera.viewport_to_world(camera_transform, cursor))
        .map(|ray| ray.origin.truncate())
    {
        let world_position_vec3 = Vec3::new(world_position.x, world_position.y, 0.0);
        mycoords.0 = world_position_vec3;
        return Some(world_position_vec3);
    }

    // if the cursor is not inside the window, we don't update the world position
    return None;
}

fn get_primary_window_size(q_window: &Query<&Window, With<PrimaryWindow>>) -> Vec2 {
    let window = q_window.get_single().unwrap();
    let window = Vec2::new(window.width() as f32, window.height() as f32);
    window
}