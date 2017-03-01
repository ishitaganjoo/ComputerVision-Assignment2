#ifndef __RANDPOINTS_H__
#define __RANDPOINTS_H__


class randPoints{



public:
    static float _x;
    static float _y;
    static float _xP;
    static float _yP;
randPoints() {
  }

  randPoints(float x, float y, float xP, float yP){
    _x = x;
    _y = y;
    _xP  = xP;
    _yP = yP;
  }

};

#endif
