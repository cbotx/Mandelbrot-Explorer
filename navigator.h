#ifndef __NAVIGATOR_H__
#define __NAVIGATOR_H__

#include <chrono>
#include <algorithm>
#include <mutex>
#include <array>

class Navigator {
protected:
    int _w, _h;
    double _dx, _dy;
    double _display_dx, _display_dy;
    double _k;

private:
    double _zoom_progress;
    double _zoom_start, _zoom_end;
    double _zoom_step;
    double _zoom_time; // zoom animation time in ms
    double _zoom_x, _zoom_y; // zoom center
    
    std::mutex _mu;
    int _drag_old_x, _drag_old_y, _drag_unprocessed_dx, _drag_unprocessed_dy;
    bool _dragging;
    
    std::chrono::time_point<std::chrono::high_resolution_clock> _start_time;

    void (*_fix_image_cb)(void);

public:
    Navigator(int width, int height, double zoom_step, double zoom_time);

    virtual ~Navigator() = default;

    double easeFunction(double l, double r, double t);

    void ZoomIn(int x, int y);

    void ZoomOut(int x, int y);

    void DragStart(int x, int y);

    void Drag(int x, int y);

    void DragEnd();

    void Update();

    void BindFixImageCallback(void (*func)(void));

    std::array<std::array<double, 2>, 4> GetBoundary();

private:
    void init();
};

#endif