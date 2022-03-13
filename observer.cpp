#include "observer.h"

void ColorSubject::Attach(ColorObserver* observer) {
    _observers.push_back(observer);
}

void ColorSubject::Notify(float R, float G, float B) {
    for (auto& observer : _observers) {
        observer->UpdateColor(R, G, B);
    }
}
