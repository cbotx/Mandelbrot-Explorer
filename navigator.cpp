#include <iostream>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <math.h>

#include "navigator.h"

Navigator::Navigator(int width, int height, double zoom_step, double zoom_time):
        _w(width), _h(height), _zoom_step(zoom_step), _zoom_time(zoom_time) {
    _fix_image_cb = nullptr;
    init();
}

double Navigator::easeFunction(double l, double r, double t) {
    return (r - l) * std::min(std::max(t, 0.0), 1.0) + l;
}

void Navigator::ZoomIn(int x, int y) {
    y = _h - y;
    if (_zoom_progress < 1.0) {
        _dx += (x - _zoom_x) / _k * (_k - 1);
        _dy += (y - _zoom_y) / _k * (_k - 1);
    }
    _zoom_x = x;
    _zoom_y = y;
    _zoom_start = easeFunction(_zoom_start, _zoom_end, _zoom_progress);
    _zoom_progress = 0.0;
    _zoom_end += _zoom_step;
    _start_time = std::chrono::high_resolution_clock::now();
}

void Navigator::ZoomOut(int x, int y) {
    y = _h - y;
    if (_zoom_progress < 1.0) {
        _dx += (x - _zoom_x) / _k * (_k - 1);
        _dy += (y - _zoom_y) / _k * (_k - 1);
    }
    _zoom_x = x;
    _zoom_y = y;
    _zoom_start = easeFunction(_zoom_start, _zoom_end, _zoom_progress);
    _zoom_progress = 0.0;
    _zoom_end -= _zoom_step;
    _start_time = std::chrono::high_resolution_clock::now();
}

void Navigator::DragStart(int x, int y) {
    y = _h - y;
    _drag_old_x = x;
    _drag_old_y = y;
    _dragging = true;
}

void Navigator::Drag(int x, int y) {
    y = _h - y;
    {
        const std::lock_guard<std::mutex> lock(_mu);
        _drag_unprocessed_dx += x - _drag_old_x;
        _drag_unprocessed_dy += y - _drag_old_y;
    }
    _drag_old_x = x;
    _drag_old_y = y;
}

void Navigator::DragEnd() {
    _dragging = false;
    if (_zoom_progress >= 1.0) {
        init();
    }
}

void Navigator::Update() {
    if (_zoom_progress < 1.0) {
        std::chrono::time_point<std::chrono::high_resolution_clock> end_time = std::chrono::high_resolution_clock::now();
        int elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - _start_time).count();
        _zoom_progress = 1.0 * elapsed_time / _zoom_time;
        if (_zoom_progress >= 1.0) {
            _zoom_progress = 1.0;
            if (!_dragging) {
                init();
            }
        }
    }
    _k = std::pow(2.0, easeFunction(_zoom_start, _zoom_end, _zoom_progress));
    {
        const std::lock_guard<std::mutex> lock(_mu);
        _dx += _drag_unprocessed_dx / _k;
        _dy += _drag_unprocessed_dy / _k;
        _drag_unprocessed_dx = 0;
        _drag_unprocessed_dy = 0;
    }
    _display_dx = _dx - (_zoom_x - _dx) * (_k - 1);
    _display_dy = _dy - (_zoom_y - _dy) * (_k - 1);

}
void Navigator::init() {
    if (_fix_image_cb) _fix_image_cb();
    _dx = _dy = _display_dx = _display_dy = 0;
    _drag_unprocessed_dx = _drag_unprocessed_dy = 0;
    _zoom_start = _zoom_end = 0;
    _zoom_x = _zoom_y = 0;
    _k = 1;
    _zoom_progress = 1;
    _dragging = false;
}


void Navigator::BindFixImageCallback(void (*func)(void)) {
    _fix_image_cb = func;
}

std::array<std::array<double, 2>, 4> Navigator::GetBoundary() {
    double lx = -_display_dx / (_k * _w);
    double ly = -_display_dy / (_k * _h);
    double rx = (_w - _display_dx) / (_k * _w);
    double ry = (_h - _display_dy) / (_k * _h);
    std::array<std::array<double, 2>, 4> arr = {lx, ly, lx, ry, rx, ry, rx, ly};
    return arr;
}