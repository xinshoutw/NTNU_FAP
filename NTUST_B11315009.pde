/**
 * @file NTUST_B11315009.pde
 * @author xinshoutw <contact@xinshou.tw>
 * @date 2025/03/10
 * @version 0.0.1
 * 
 * @brief Animated art sketch with circles and lines transition effects.
 *
 * This program combines the entry, exit, and animation effects of circles and lines,
 * using Doxygen format to provide detailed explanations of constants and functions.
 */

// ===== Global Animation Parameters & Constants =====

/**
 * @brief Duration of circle entry (falling from top) in milliseconds
 */
final int CIRCLE_ENTER_DURATION = 800;

/**
 * @brief Total duration of circle scaling (expansion/contraction) effect in milliseconds
 */
final int CIRCLE_SCALE_DURATION = 500;

/**
 * @brief Duration of circle exit (moving downward off-screen) in milliseconds
 */
final int CIRCLE_EXIT_DURATION = 300;

/**
 * @brief Duration of drawing each line in milliseconds
 */
final int LINE_DRAW_DURATION = 200;

/**
 * @brief Duration of line exit (gradual disappearance from starting point) in milliseconds
 */
final int LINE_EXIT_DURATION = 300;

/**
 * @brief Rotation speed of circles after settling, in radians per millisecond
 */
final float CIRCLE_ROTATION_SPEED = 0.0001 * TWO_PI;

/**
 * @brief Maximum scaling factor of circles after expansion
 */
final float CIRCLE_SCALE_MAX = 1.4;

/**
 * @brief Maximum scaling factor of circles after expansion
 */
final float CIRCLE_BOUNCE_OFFSET = 0.75;

ArrayList < CirclePiece > circlePieces; ///< Stores circle/arc segments
ArrayList < LineObj > lines; ///< Stores various lines

int artStartTime; ///< Start time of entry animation in milliseconds
boolean isExiting = false; ///< Whether the exit animation is currently active
int exitStartTime; ///< Start time of exit animation in milliseconds

// ------------------------------------------------------

/**
 * @brief Processing setup function
 */
void setup() {
  size(800, 600);
  generateArt();
}

/**
 * @brief Main drawing loop
 *
 * Calls either the normal animation or exit animation function based on the exit state
 */
void draw() {
  drawGradientBG(color(230, 245, 255), color(180, 210, 240));
  if (isExiting) {
    drawExitAnimation();
  } else {
    drawNormalAnimation();
  }
}

/**
 * @brief Mouse click event handler
 *
 * Starts the exit animation if not already in the exit state
 */
void mousePressed() {
  if (!isExiting) {
    isExiting = true;
    exitStartTime = millis();
  }
  saveFrame();
}

/**
 * @brief Generates new artwork: initializes circles and lines and resets animation start time
 */
void generateArt() {
  circlePieces = new ArrayList < CirclePiece > ();
  lines = new ArrayList < LineObj > ();
  generateLines();
  generateCircles();
  artStartTime = millis();
}

/**
 * @brief Randomly generates lines and adds them to the global lines list
 */
void generateLines() {
  int lineCount = int(random(4, 8));
  for (int i = 0; i < lineCount; i++) {
    float t = random(1);
    if (t < 0.5) {
      lines.add(new StraightLine());
    } else if (t < 0.8) {
      lines.add(new CurvedLine());
    } else {
      lines.add(new ZigZagLine());
    }
  }
}

/**
 * @brief Randomly generates circles, splits them with lines, and adds them to the global circle list
 */
void generateCircles() {
  int circleNum = int(random(10, 16));
  int attempts = 0;
  while (circlePieces.size() < circleNum && attempts < circleNum * 50) {
    float r = random(60, 130);
    float x = random(r, width - r);
    float y = random(r, height - r);
    if (!tooOverlapping(x, y, r)) {
      color c = color(random(70, 120), random(130, 180), random(200, 255));
      CirclePiece base = new CirclePiece(x, y, r, 0, TWO_PI, c);
      ArrayList < CirclePiece > splitted = splitCircleWithLines(base, lines);
      circlePieces.addAll(splitted);
    }
    attempts++;
  }
}

/**
 * @brief Checks if a new circle overlaps excessively with existing circles
 * @param x X-coordinate of the new circle's center
 * @param y Y-coordinate of the new circle's center
 * @param r Radius of the new circle
 * @return true if excessively overlapping, false otherwise
 */
boolean tooOverlapping(float x, float y, float r) {
  for (CirclePiece cp: circlePieces) {
    if (dist(x, y, cp.cx, cp.cy) < (r + cp.r) * 0.6) {
      return true;
    }
  }
  return false;
}

/**
 * @brief Splits a circle with lines and returns the resulting arc segments
 * @param base The original circle
 * @param lines The collection of lines used for splitting
 * @return The collection of resulting arc segments
 */
ArrayList < CirclePiece > splitCircleWithLines(CirclePiece base, ArrayList < LineObj > lines) {
  ArrayList < CirclePiece > arcs = new ArrayList < CirclePiece > ();
  arcs.add(base);
  for (LineObj ln: lines) {
    ArrayList < CirclePiece > newArcs = new ArrayList < CirclePiece > ();
    for (CirclePiece arc: arcs) {
      float[] interAngles = ln.getCircleIntersections(arc.cx, arc.cy, arc.r);
      if (interAngles == null || interAngles.length == 0) {
        newArcs.add(arc);
        continue;
      }
      ArrayList < Float > cutAngles = new ArrayList < Float > ();
      for (float angle: interAngles) {
        angle = (angle % TWO_PI + TWO_PI) % TWO_PI;
        if (angleInRange(angle, arc.startAng, arc.endAng)) {
          cutAngles.add(angle);
        }
      }
      if (cutAngles.size() == 0) {
        newArcs.add(arc);
      } else {
        cutAngles.sort(Float::compareTo);
        float prev = arc.startAng;
        for (float cAng: cutAngles) {
          if (cAng > prev) {
            newArcs.add(arc.makeSubArc(prev, cAng));
            prev = cAng;
          }
        }
        if (prev < arc.endAng) {
          newArcs.add(arc.makeSubArc(prev, arc.endAng));
        }
      }
    }
    arcs = newArcs;
  }
  for (CirclePiece cp: arcs) {
    if ((cp.endAng - cp.startAng) < TWO_PI) {
      float delta = random(-20, 20);
      cp.c = color(red(cp.c) + delta, green(cp.c) + delta, blue(cp.c) + delta, 255);
    }
  }
  return arcs;
}

/**
 * @brief Determines if an angle falls within a specified range
 * @param angle The target angle
 * @param start The start angle of the range
 * @param end The end angle of the range
 * @return true if the angle is within [start, end], false otherwise
 */
boolean angleInRange(float angle, float start, float end) {
  return (angle >= start && angle <= end);
}

/**
 * @brief Draws a gradient background effect
 * @param c1 Starting color of the gradient
 * @param c2 Ending color of the gradient
 */
void drawGradientBG(color c1, color c2) {
  for (int y = 0; y < height; y++) {
    stroke(lerpColor(c1, c2, map(y, 0, height, 0, 1)));
    line(0, y, width, y);
  }
  noStroke();
}

/**
 * @brief Cubic ease-out easing function
 * @param t Progress (0 to 1)
 * @return Eased value
 */
float easeOut(float t) {
  return 1 - pow(1 - t, 3);
}

/**
 * @brief Calculates circle scaling factor based on progress
 * @param t Progress (0 to 1)
 * @return Scaling factor
 */
float getCircleScale(float t) {
  return 1 + (CIRCLE_SCALE_MAX - 1) * sin(PI * t);
}

/**
 * @brief Gets the scaling factor upon landing, using a damped oscillation function to produce multiple natural oscillations.
 * 
 * This function uses an exponentially damped oscillation formula, with an angular frequency set to (e * TWO_PI),
 * linking the oscillation frequency to the natural number e, without a fixed number of oscillations,
 * until the scaling effect gradually decays.
 * 
 * @param t Landing scaling animation progress (0 to 1)
 * @return Current scaling factor
 */
float getLandingScale(float t) {
  float A = CIRCLE_SCALE_MAX - 1.0; // Initial amplitude
  float lambda = 3.0; // Damping coefficient, adjustable to control oscillation decay speed
  float omega = (float) Math.E * TWO_PI; // Angular frequency, related to the natural number e
  return 1.0 + A * exp(-lambda * t) * cos(omega * t);
}

/**
 * @brief Normal animation: circle entry and line-by-line drawing
 */
void drawNormalAnimation() {
  int elapsed = millis() - artStartTime;
  float circleEnterProg = constrain(elapsed / (float) CIRCLE_ENTER_DURATION, 0, 1);
  float rotationAngle = 0;
  if (circleEnterProg >= 1) {
    rotationAngle = CIRCLE_ROTATION_SPEED * (millis() - (artStartTime + CIRCLE_ENTER_DURATION));
  }

  // Draw circle shadows and bodies
  for (CirclePiece cp: circlePieces) {
    cp.displayShadowEnter(circleEnterProg, rotationAngle);
    cp.displayEnter(circleEnterProg, rotationAngle);
  }

  // Draw lines (after circle entry animation completes)
  if (elapsed >= CIRCLE_ENTER_DURATION) {
    int lineElapsed = millis() - (artStartTime + CIRCLE_ENTER_DURATION);
    int currentLineIndex = lineElapsed / LINE_DRAW_DURATION;
    float lineDrawProg = (lineElapsed % LINE_DRAW_DURATION) / (float) LINE_DRAW_DURATION;
    for (int i = 0; i < lines.size(); i++) {
      if (i < currentLineIndex) {
        lines.get(i).display();
      } else if (i == currentLineIndex) {
        lines.get(i).displayAnimated(lineDrawProg);
      }
    }
  }
}

/**
 * @brief Exit animation: circles move downward and lines gradually disappear from the start
 */
void drawExitAnimation() {
  float circleExitProg = constrain((millis() - exitStartTime) / (float) CIRCLE_EXIT_DURATION, 0, 1);
  float lineExitProg = constrain((millis() - exitStartTime) / (float) LINE_EXIT_DURATION, 0, 1);

  for (CirclePiece cp: circlePieces) {
    cp.displayShadowExit(circleExitProg);
    cp.displayExit(circleExitProg);
  }
  for (LineObj ln: lines) {
    ln.displayExit(lineExitProg);
  }

  if (circleExitProg >= 1 && lineExitProg >= 1) {
    generateArt();
    isExiting = false;
  }
}

// ==========================================================
// ===================== CirclePiece Class ==================
// ==========================================================

/**
 * @brief Represents a circular arc (circle segment)
 */
class CirclePiece {
  float cx, cy; ///< Center coordinates
  float r; ///< Radius
  float startAng, endAng; ///< Start and end angles of the arc
  color c; ///< Color

  /**
   * @brief Constructor
   * @param cx X-coordinate of the center
   * @param cy Y-coordinate of the center
   * @param r Radius
   * @param startAng Start angle (radians)
   * @param endAng End angle (radians)
   * @param c Color
   */
  CirclePiece(float cx, float cy, float r, float startAng, float endAng, color c) {
    this.cx = cx;
    this.cy = cy;
    this.r = r;
    this.startAng = min(startAng, endAng);
    this.endAng = max(startAng, endAng);
    this.c = c;
  }

  /**
   * @brief Creates a sub-arc from the current arc
   * @param s New start angle
   * @param e New end angle
   * @return A new CirclePiece object
   */
  CirclePiece makeSubArc(float s, float e) {
    return new CirclePiece(cx, cy, r, s, e, c);
  }

  /**
   * @brief Entry animation: draws the shadow (including fall and scaling effects)
   * @param prog Animation progress (0 to 1)
   * @param rot Rotation angle (radians)
   */
  void displayShadowEnter(float prog, float rot) {
    pushMatrix();
    if (prog < 1) {
      float currentY = lerp(-r * 2, cy, easeOut(prog));
      float scaleProg = min(prog * CIRCLE_ENTER_DURATION / (float) CIRCLE_SCALE_DURATION, 1);
      float s = getCircleScale(scaleProg);
      translate(cx + 5, currentY + 5);
      scale(s);
    } else {
      translate(cx + 5, cy + 5);
      rotate(rot);
    }
    fill(0, 50);
    noStroke();
    arc(0, 0, r * 2, r * 2, startAng, endAng);
    popMatrix();
  }

  /**
   * @brief Modified circle entry animation: triggers landing scaling oscillation effect after falling
   * 
   * When animation progress prog is less than 1, the circle only performs the falling animation;
   * when prog >= 1, the circle stays in its final position and begins the landing scaling oscillation animation,
   * with the scaling effect determined by the getLandingScale() function, combined with rotation.
   * 
   * @param prog Animation progress (0 to 1)
   * @param rot Rotation angle (radians)
   */
  void displayEnter(float prog, float rot) {
    pushMatrix();
    float effectiveProg = min(prog, 1);
    float currentY = lerp(-r * 2, cy, easeOut(effectiveProg));
    translate(cx, currentY);

    float blend = constrain((prog - CIRCLE_BOUNCE_OFFSET) / 0.1, 0, 1);
    if (blend > 0) {
      float landingTime = millis() - (artStartTime + CIRCLE_ENTER_DURATION);
      float landingProg = constrain(landingTime / (float) CIRCLE_SCALE_DURATION, 0, 1);
      float s = getLandingScale(landingProg);

      scale(lerp(1, s, blend));
      rotate(lerp(0, rot, blend));
    }

    fill(c);
    noStroke();
    arc(0, 0, r * 2, r * 2, startAng, endAng);
    popMatrix();
  }

  /**
   * @brief Exit animation: draws the shadow (only downward movement)
   * @param prog Exit progress (0 to 1)
   */
  void displayShadowExit(float prog) {
    pushMatrix();
    translate(cx + 5, cy + 5 + lerp(0, height, prog));
    fill(0, 50);
    noStroke();
    arc(0, 0, r * 2, r * 2, startAng, endAng);
    popMatrix();
  }

  /**
   * @brief Exit animation: draws the circle body (moves downward off-screen)
   * @param prog Exit progress (0 to 1)
   */
  void displayExit(float prog) {
    pushMatrix();
    translate(cx, cy + lerp(0, height, prog));
    fill(c);
    noStroke();
    arc(0, 0, r * 2, r * 2, startAng, endAng);
    popMatrix();
  }
}

// ==========================================================
// ==================== LineObj Base Class ==================
// ==========================================================

/**
 * @brief Abstract base class for line objects
 */
abstract class LineObj {
  color lineColor; ///< Line color
  float sw; ///< Stroke width

  /**
   * @brief Constructor, randomly sets color and stroke width
   */
  LineObj() {
    lineColor = color(random(30, 80), random(60, 160), random(180, 255));
    sw = random(1, 3);
  }

  /**
   * @brief Calculates intersections between the line and a specified circle (in polar angles on the circle)
   * @param cx X-coordinate of the circle's center
   * @param cy Y-coordinate of the circle's center
   * @param r Circle radius
   * @return Array of intersection angles, or null if no intersections
   */
  abstract float[] getCircleIntersections(float cx, float cy, float r);

  /**
   * @brief Draws the complete line
   */
  abstract void display();

  /**
   * @brief Draws the line animation (partial drawing)
   * @param prog Animation progress (0 to 1)
   */
  abstract void displayAnimated(float prog);

  /**
   * @brief Draws the line exit animation (gradual disappearance from the start)
   * @param prog Exit progress (0 to 1)
   */
  abstract void displayExit(float prog);
}

// ------------------------------------------------------
// ================= StraightLine Class =================
// ------------------------------------------------------

/**
 * @brief Straight line object
 */
class StraightLine extends LineObj {
  float x1, y1, x2, y2; ///< Coordinates of the line's endpoints

  /**
   * @brief Constructor, randomly generates line endpoints
   */
  StraightLine() {
    super();
    x1 = random(width);
    y1 = random(height);
    x2 = random(width);
    y2 = random(height);
  }

  @Override
  void display() {
    stroke(lineColor);
    strokeWeight(sw);
    line(x1, y1, x2, y2);
    noStroke();
  }

  @Override
  void displayAnimated(float prog) {
    stroke(lineColor);
    strokeWeight(sw);
    float x = lerp(x1, x2, prog);
    float y = lerp(y1, y2, prog);
    line(x1, y1, x, y);
    noStroke();
  }

  @Override
  void displayExit(float prog) {
    stroke(lineColor);
    strokeWeight(sw);
    float newX1 = lerp(x1, x2, prog);
    float newY1 = lerp(y1, y2, prog);
    line(newX1, newY1, x2, y2);
    noStroke();
  }

  @Override
  float[] getCircleIntersections(float cx, float cy, float r) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float fx = x1 - cx;
    float fy = y1 - cy;
    float A = dx * dx + dy * dy;
    float B = 2 * (fx * dx + fy * dy);
    float C = fx * fx + fy * fy - r * r;
    float disc = B * B - 4 * A * C;
    if (disc < 0) return null;
    ArrayList < Float > angles = new ArrayList < Float > ();
    if (disc == 0) {
      float t = -B / (2 * A);
      float ix = x1 + t * dx;
      float iy = y1 + t * dy;
      angles.add(atan2(iy - cy, ix - cx));
    } else {
      float sqrtDisc = sqrt(disc);
      float t1 = (-B + sqrtDisc) / (2 * A);
      float t2 = (-B - sqrtDisc) / (2 * A);
      float ix1 = x1 + t1 * dx;
      float iy1 = y1 + t1 * dy;
      float ix2 = x1 + t2 * dx;
      float iy2 = y1 + t2 * dy;
      angles.add(atan2(iy1 - cy, ix1 - cx));
      angles.add(atan2(iy2 - cy, ix2 - cx));
    }
    float[] res = new float[angles.size()];
    for (int i = 0; i < angles.size(); i++) {
      res[i] = angles.get(i);
    }
    return res;
  }
}

// ------------------------------------------------------
// ================= CurvedLine Class =================
// ------------------------------------------------------

/**
 * @brief Curved line object, implemented with a quadratic Bézier curve
 */
class CurvedLine extends LineObj {
  PVector p1, p2, p3; ///< Control points of the curve

  /**
   * @brief Constructor, randomly generates control points
   */
  CurvedLine() {
    super();
    p1 = new PVector(random(width), random(height));
    p2 = new PVector(random(width), random(height));
    p3 = new PVector(random(width), random(height));
  }

  @Override
  void display() {
    noFill();
    stroke(lineColor);
    strokeWeight(sw);
    beginShape();
    vertex(p1.x, p1.y);
    quadraticVertex(p2.x, p2.y, p3.x, p3.y);
    endShape();
    noStroke();
  }

  @Override
  void displayAnimated(float prog) {
    noFill();
    stroke(lineColor);
    strokeWeight(sw);
    beginShape();
    int steps = max(2, int(20 * prog));
    for (int i = 0; i <= steps; i++) {
      float t = i / (float) steps;
      if (t > prog) break;
      PVector p = getQuadPoint(t);
      vertex(p.x, p.y);
    }
    endShape();
    noStroke();
  }

  @Override
  void displayExit(float prog) {
    noFill();
    stroke(lineColor);
    strokeWeight(sw);
    beginShape();
    int steps = 20;
    for (int i = 0; i <= steps; i++) {
      float t = lerp(prog, 1, i / (float) steps);
      PVector p = getQuadPoint(t);
      vertex(p.x, p.y);
    }
    endShape();
    noStroke();
  }

  /**
   * @brief Gets the point on the Bézier curve corresponding to parameter t
   * @param t Parameter (0 to 1)
   * @return Point on the curve
   */
  PVector getQuadPoint(float t) {
    float u = 1 - t;
    float x = u * u * p1.x + 2 * u * t * p2.x + t * t * p3.x;
    float y = u * u * p1.y + 2 * u * t * p2.y + t * t * p3.y;
    return new PVector(x, y);
  }

  @Override
  float[] getCircleIntersections(float cx, float cy, float r) {
    int SAMPLES = 200;
    ArrayList < Float > angles = new ArrayList < Float > ();
    PVector prev = getQuadPoint(0);
    for (int i = 1; i <= SAMPLES; i++) {
      float t = i / (float) SAMPLES;
      PVector curr = getQuadPoint(t);
      float[] subAngles = lineCircleInter(prev.x, prev.y, curr.x, curr.y, cx, cy, r);
      if (subAngles != null) {
        for (float a: subAngles) angles.add(a);
      }
      prev = curr;
    }
    if (angles.size() == 0) return null;
    float[] res = new float[angles.size()];
    for (int i = 0; i < angles.size(); i++) {
      res[i] = angles.get(i);
    }
    return res;
  }
}

// ------------------------------------------------------
// ================= ZigZagLine Class =================
// ------------------------------------------------------

/**
 * @brief Zigzag line object, using multiple points to create a jagged effect
 */
class ZigZagLine extends LineObj {
  ArrayList < PVector > pts; ///< Collection of points forming the zigzag line

  /**
   * @brief Constructor, randomly generates 3 to 5 points
   */
  ZigZagLine() {
    super();
    int n = int(random(3, 6));
    pts = new ArrayList < PVector > ();
    for (int i = 0; i < n; i++) {
      pts.add(new PVector(random(width), random(height)));
    }
  }

  @Override
  void display() {
    stroke(lineColor);
    strokeWeight(sw);
    noFill();
    beginShape();
    for (PVector v: pts) {
      vertex(v.x, v.y);
    }
    endShape();
    noStroke();
  }

  @Override
  void displayAnimated(float prog) {
    stroke(lineColor);
    strokeWeight(sw);
    noFill();
    beginShape();
    int totalSegments = pts.size() - 1;
    float totalProgress = prog * totalSegments;
    int segFullyDrawn = floor(totalProgress);
    for (int i = 0; i <= segFullyDrawn && i < pts.size(); i++) {
      vertex(pts.get(i).x, pts.get(i).y);
    }
    if (segFullyDrawn < totalSegments) {
      float segProg = totalProgress - segFullyDrawn;
      PVector a = pts.get(segFullyDrawn);
      PVector b = pts.get(segFullyDrawn + 1);
      vertex(lerp(a.x, b.x, segProg), lerp(a.y, b.y, segProg));
    }
    endShape();
    noStroke();
  }

  @Override
  void displayExit(float prog) {
    stroke(lineColor);
    strokeWeight(sw);
    noFill();
    beginShape();
    int totalSegments = pts.size() - 1;
    float totalProgress = prog * totalSegments;
    int segIndex = floor(totalProgress);
    float segFrac = totalProgress - segIndex;
    if (segIndex < pts.size() - 1) {
      PVector a = pts.get(segIndex);
      PVector b = pts.get(segIndex + 1);
      vertex(lerp(a.x, b.x, segFrac), lerp(a.y, b.y, segFrac));
      for (int i = segIndex + 1; i < pts.size(); i++) {
        vertex(pts.get(i).x, pts.get(i).y);
      }
    }
    endShape();
    noStroke();
  }

  @Override
  float[] getCircleIntersections(float cx, float cy, float r) {
    ArrayList < Float > angles = new ArrayList < Float > ();
    for (int i = 0; i < pts.size() - 1; i++) {
      PVector a = pts.get(i);
      PVector b = pts.get(i + 1);
      float[] subAngles = lineCircleInter(a.x, a.y, b.x, b.y, cx, cy, r);
      if (subAngles != null) {
        for (float aa: subAngles) angles.add(aa);
      }
    }
    if (angles.size() == 0) return null;
    float[] res = new float[angles.size()];
    for (int i = 0; i < angles.size(); i++) {
      res[i] = angles.get(i);
    }
    return res;
  }
}

// ------------------------------------------------------
// ================= Utility Functions ==================
// ------------------------------------------------------

/**
 * @brief Calculates intersections between a line segment and a circle, expressed as polar angles on the circle
 * @param x1 X-coordinate of the line segment's start point
 * @param y1 Y-coordinate of the line segment's start point
 * @param x2 X-coordinate of the line segment's end point
 * @param y2 Y-coordinate of the line segment's end point
 * @param cx X-coordinate of the circle's center
 * @param cy Y-coordinate of the circle's center
 * @param r Circle radius
 * @return Array of intersection angles (radians), or null if no intersections
 */
float[] lineCircleInter(float x1, float y1, float x2, float y2, float cx, float cy, float r) {
  float dx = x2 - x1;
  float dy = y2 - y1;
  float fx = x1 - cx;
  float fy = y1 - cy;
  float A = dx * dx + dy * dy;
  float B = 2 * (fx * dx + fy * dy);
  float C = fx * fx + fy * fy - r * r;
  float disc = B * B - 4 * A * C;
  if (disc < 0) return null;
  ArrayList < Float > angles = new ArrayList < Float > ();
  if (disc == 0) {
    float t = -B / (2 * A);
    if (t >= 0 && t <= 1) {
      float ix = x1 + t * dx;
      float iy = y1 + t * dy;
      angles.add(atan2(iy - cy, ix - cx));
    }
  } else {
    float sqrtDisc = sqrt(disc);
    float t1 = (-B + sqrtDisc) / (2 * A);
    float t2 = (-B - sqrtDisc) / (2 * A);
    if (t1 >= 0 && t1 <= 1) {
      float ix1 = x1 + t1 * dx;
      float iy1 = y1 + t1 * dy;
      angles.add(atan2(iy1 - cy, ix1 - cx));
    }
    if (t2 >= 0 && t2 <= 1) {
      float ix2 = x1 + t2 * dx;
      float iy2 = y1 + t2 * dy;
      angles.add(atan2(iy2 - cy, ix2 - cx));
    }
  }
  if (angles.size() == 0) return null;
  float[] res = new float[angles.size()];
  for (int i = 0; i < angles.size(); i++) {
    res[i] = angles.get(i);
  }
  return res;
}
