#include <stdint.h>
#include <iostream>
#include <mpir.h>
#include <thread>
#include <mutex>
#include <chrono>

#include "glut_canvas.h"
#include "navigator.h"
#include "mandel_perturbation.h"
#include "mandel_navigator.h"
#include "Image.h"
#include "gradient_window.h"

constexpr int WIDTH = 600;
constexpr int HEIGHT = 400;
constexpr int FPS = 60;
constexpr int MAX_ITERATION = 1000000;

GlutCanvas* canvas;
GlutCanvas::Texture texture;
uint8_t* bitmap;
MandelNavigator* nav;
bool mouse_down;
float* iter;
std::mutex mu;

void mouse_move(int x, int y) {
    if (mouse_down) {
        nav->Drag(x, y);
    }
}

void mouse_click(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        if (button == GLUT_LEFT_BUTTON) {
            mouse_down = true;
            nav->DragStart(x, y);
        }
    } else if (state == GLUT_UP) {
        if (button == GLUT_LEFT_BUTTON) {
            mouse_down = false;
            nav->DragEnd();
        } else if (button == 3) {
            nav->ZoomIn(x, y);
        } else if (button == 4) {
            nav->ZoomOut(x, y);
        }
    }
}


void key_press(unsigned char key, int x, int y) {
    if (key == 'r' || key == 'R') {
        nav->Reset();
    }
}

void display() {
    static std::chrono::time_point<std::chrono::high_resolution_clock> last_time = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> cur_time = std::chrono::high_resolution_clock::now();
    int dt = std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - last_time).count();
    if (dt > 400) {
        last_time = cur_time;
        nav->UpdateBitmap(bitmap);
        canvas->BitmapToTexture(texture, bitmap);
    }
    
    nav->Update();
    canvas->BeginDraw();
    canvas->RenderTexture(texture, nav->GetBoundary(), {-1, -1, -1, 1, 1, 1, 1, -1});
    canvas->EndDraw();
}

void reshape(int width, int height) {
    canvas->GlutResetWindowSize();
}

void fix_image() {
    {
        const std::lock_guard<std::mutex> lock(mu);
        canvas->ScreenToBitmap(bitmap);
        canvas->BitmapToTexture(texture, bitmap);
    }
    nav->UpdateCoords();
    nav->StartCompute();
}

// void refresh_loop(int interval) {
//     while(true) {
//         std::this_thread::sleep_for(std::chrono::milliseconds(interval));
//         {
//             const std::lock_guard<std::mutex> lock(mu);
//             nav->UpdateBitmap(bitmap);
//             if (texture) canvas->BitmapToTexture(texture, bitmap);
//         }
//     }
// }

int main(int argc, char** argv) {
    bitmap = new uint8_t[WIDTH * HEIGHT * 3];

    nav = new MandelNavigator(WIDTH, HEIGHT, 3, 1000000, 1, 500);
    nav->BindFixImageCallback(fix_image);

    GlutCanvas::GlutInit(argc, argv);
    canvas = new GlutCanvas();
    int mandel_window = canvas->GlutCreateWindow(WIDTH, HEIGHT, 100, 100, "Mandel  |  R = Reset");
    texture = canvas->CreateTexture();
    canvas->GlutBindMotionFunc(mouse_move);
    canvas->GlutBindMouseFunc(mouse_click);
    canvas->GlutBindDisplayFunc(display);
    canvas->GlutBindReshapeFunc(reshape);
    canvas->GlutBindKeyboardFunc(key_press);
    canvas->GlutSetRenderInterval(1000 / FPS);
    canvas->BitmapToTexture(texture, bitmap);

    gradient_window = new GradientWindow(770, 300, "Control Panel");
    gradient_window->BindNavigator(nav);
    gradient_window->_windows.push_back(mandel_window);

    nav->StartCompute();

    canvas->GlutRun();
    return 0;
}