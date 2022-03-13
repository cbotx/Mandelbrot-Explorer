#ifndef __OBSERVEr_H__
#define __OBSERVER_H__

#include <vector>

class ColorObserver {

    friend class ColorSubject;

protected:
    virtual void UpdateColor(float R, float G, float B);
};

class ColorSubject {
public:
  void Attach(ColorObserver* observer);

  void Notify(float R, float G, float B);

private:
    std::vector<ColorObserver* > _observers;
};

#endif