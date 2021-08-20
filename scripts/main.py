import cv2
import numpy as np
import networkx
import click
import math


@click.command()
@click.option(
    "--image-path",
    required=True,
    default="",
    help="Specify the path to the bondgraph image"
)
def main(image_path):
    img = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    # Perform Canny Transform
    dst = cv2.Canny(img, 50, 200, None, 3)
    # Copy edges to the images that will display the results in the BGR
    cdst = cv2.cvtColor(dst, cv2.COLOR_GRAY2BGR)
    cdstP = np.copy(cdst)
    # Find the lines
    lines = cv2.HoughLines(dst, 1, np.pi/100, 150, None, 0, 0)

    if lines is not None:
        for i in range(0, len(lines)):
            rho = lines[i][0][0]
            theta = lines[i][0][1]
            a = np.cos(theta)
            b = np.sin(theta)
            x0 = a * rho
            y0 = b * rho
            pt1 = (int(x0 + 1000*(-b)), int(y0 + 1000*(a)))
            pt2 = (int(x0 - 1000*(-b)), int(y0 - 1000*(a)))
            cv2.line(cdst, pt1, pt2, (0, 0, 255), 3, cv2.LINE_AA)
    linesP = cv2.HoughLinesP(dst, 1, 0.1*np.pi / 180, 50, None, 50, 10)
    if linesP is not None:
        for i in range(0, len(linesP)):
            l = linesP[i][0]
            cv2.line(cdstP, (l[0], l[1]), (l[2], l[3]),
                     (0, 0, 255), 3, cv2.LINE_AA)
    cv2.imshow("Source", cv2.resize(img, (0,0), fx=0.2, fy=0.2))
    # cv2.imshow("Detected Lines (in red) - Standard Hough Line Transform",
            #    cv2.resize(cdst, (0, 0), fx=0.2, fy=0.2))
    cv2.imshow("Detected Lines (in red) - Probabilistic Line Transform",
               cv2.resize(cdstP, (0, 0), fx=0.2, fy=0.2))
    if cv2.waitKey(0) == ord('q'):
        _
    return 0


if __name__ == "__main__":
    main()
