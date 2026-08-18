#ifndef GUI_EXPORT_H
#define GUI_EXPORT_H

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <cstdint>
#include <string>
#include <vector>

class MandelNavigator;

// Open the modal high-resolution PNG export dialog for the current view.
void showExportDialog(HWND owner, MandelNavigator* nav);
// Helpers reused by the main window's quick-save (S key).
std::wstring defaultSaveName();
bool writeExportPNG(const wchar_t* path, const std::vector<uint8_t>& rgb, int W, int H);
bool runExportSelfTest(const wchar_t* resultPath);

#endif
