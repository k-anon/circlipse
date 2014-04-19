/** 
 * Circlipse: Ellipse-From-5-Points and Circle-From-3-Points for Illustrator
 * -------------------------------------------------------------------------
 * Given a path with five anchor points, draws an ellipse passing through all
 * given points. Given a path with three anchor points, draws a circle passing
 * through all given points.
 *
 * By K_A. Inspired by a rather different Python implementation for Inkscape
 * by epdtry.
 * v.1.0
 */
(function () {
    'use strict';

    /**
     * Given some coefficient array coeffMatrix and the matrix dimension, solves
     * while assuming the constants vector b = [1] ** n.
     */
    function solveNormalizedSystem(coeffMatrix, n) {
        var intermVector = [], solVector = [], max, pivot, i, j, k;

        /* The Doolittle Decomposition Algorithm. With it, we transform the
         * system A * x = b to L * U * x = P * b, given the lower diagonal L,
         * the upper diagonal U, and the permutation matrix P. However,
         * b = [1] ** n, so given any P derived from I_n, P * b = b. In other
         * words, no permutation of the constants vector is needed.
         */
        for (i = 0; i < n; i += 1) {
            // Apply partial pivoting here to maximize flop precision.
            max = pivot = 0;
            for (j = i; j < n; j += 1) {
                tmp = coeffMatrix[j][i];
                for (k = 0; k < i; k += 1) {
                    tmp -= coeffMatrix[j][k] * coeffMatrix[k][i];
                }
                if (Math.abs(tmp) > Math.abs(pivot)) {
                    pivot = tmp;
                    max = j;
                }
                // Optimization: Store the value into the matrix. We'll later
                // divide it by the determined pivot value (coeffMatrix[i][i]).
                coeffMatrix[j][i] = tmp;
            }
            if (!pivot) {
                // U[i][i] must not equal zero. We cannot continue.
                return null;
            }
            // Finally swap the two rows.
            if (max !== i) {
                tmp = coeffMatrix[max];
                coeffMatrix[max] = coeffMatrix[i];
                coeffMatrix[i] = tmp;
            }
            
            // Calculate lower matrix values for the bounding column by
            // dividing the predetermined dividend by the pivot value.
            for (j = i + 1; j < n; j += 1) {
                coeffMatrix[j][i] /= pivot;
            }

            // Calculate the respective row of the upper matrix.
            for (j = i + 1; j < n; j += 1) {
                tmp = coeffMatrix[i][j];
                for (k = 0; k < i; k += 1) {
                    tmp -= coeffMatrix[i][k] * coeffMatrix[k][j];
                }
                coeffMatrix[i][j] = tmp;
            }
        }

        /* We are now able to solve the equation L * U * x = P * b = b.
         * Let y = U * x to break this into two steps.
         * Step One: Solve for vector y given L * y = P * b. This can be
         * efficiently solved with the following forward-substitution algorithm.
         * L[i][i] == 1 for all integers 0 <= i < n, so this is simplified.
         */
        for (i = 0; i < n; i += 1) {
            tmp = 1.0;
            for (j = 0; j < i; j += 1) {
                tmp -= coeffMatrix[i][j] * intermVector[j];
            }
            intermVector[i] = tmp;
        }
        
        /* Step Two: Solve for the solution vector x given U * x = y. Analogous
         * to the lower-diagonal matrix L, this can be solved with the following
         * backward-substitution algorithm.
         */
        for (i = n - 1; i >= 0; i -= 1) {
            tmp = intermVector[i];
            for (j = i + 1; j < n; j += 1) {
                tmp -= coeffMatrix[i][j] * solVector[j];
            }
            solVector[i] = tmp / coeffMatrix[i][i];
        }

        return solVector;
    }

    /**
     * Given some 3-member array of two-dimensional coordinates, calculates the
     * properties of a circle passing through them. Returns null if no circle
     * is found, or an Object with the properties xOrigin, yOrigin, radius,
     * and diameter.
     */
    function calculateCircle(points) {
        var solVector, a, b, c, circleProp = {};
        
        if (points.length !== 3) {
            return null;
        }
        
        /* Consider the general conic equation:
         *
         * Ax**2 + Bxy + Cy**2 + Dx + Ey + F = 0
         *
         * For a circle, A == C, and B == 0. We can also normalize by dividing by
         * -F. Therefore, the equation of a circle can be expressed as follows
         * (with new A, B, and C):
         *
         * Ax**2 + Ay**2 + Bx + Cy = 1
         *
         * Given three points (x_i, y_i) for 0 <= i < 3, we can create a system of
         * equations to solve for all three coefficients.
         */        
        solVector = solveNormalizedSystem([
            [points[0][0] * points[0][0] + points[0][1] * points[0][1],
             points[0][0], points[0][1]],
            [points[1][0] * points[1][0] + points[1][1] * points[1][1],
             points[1][0], points[1][1]],
            [points[2][0] * points[2][0] + points[2][1] * points[2][1],
             points[2][0], points[2][1]]
        ], 3);
        
        if (!solVector) {
            return null;
        }
        
        /* Extract the coefficients from the solution vector. */
        a = solVector[0];  // x ** 2 and y ** 2 coefficient
        if (!a) {
            // That's a line, silly!
            return null;
        }
        b = solVector[1];  // x coefficient
        c = solVector[2];  // y coefficient
        
        /* Complete the squares to find (h, k)--i.e., the center point. */
        circleProp.xOrigin = -b / (2.0 * a);
        circleProp.yOrigin = -c / (2.0 * a);

        /* We complete the square on the RHS. This results in
         * 1 + A(B / (2A)) ** 2 + A(C / (2A)) ** 2
         * Normalizing the equation by dividing A gives the square of the radius
         * on the RHS.
         */
        circleProp.radius = Math.sqrt((1.0 / a + (b * b + c * c) / (4.0 * a * a)));
        if (circleProp.radius > 1e15) {
            // A rare but nasty occurrence where a is almost-but-not-quite zero,
            // resulting in a ridiculously large circle.
            return null;
        }
        circleProp.diameter = circleProp.radius * 2.0;
        
        return circleProp;
    }

    /**
     * Given some 5-member array of two-dimensional coordinates, calculates the
     * properties of an ellipse passing through them. Returns null if no
     * ellipse is found, or an Object with the properties width, height,
     * halfHeight, halfWidth, theta (rotation in radians), xOrigin, and yOrigin.
     */
    function calculateEllipse(points) {
        var solVector, a, b, c, d, e, aPrime, cPrime, dPrime, ePrime,
            theta, xOffset, yOffset, cosTheta, sinTheta, cosThetaSquared,
            sinThetaSquared, normalizer, tmp, ellipseProp = {},
            xOriginPrime, yOriginPrime, hypotenuse;

        if (points.length !== 5) {
            return null;
        }

        /* Set up a system of five linear systems based on the set of
         * five provided points. Each equation has the following form:
         *
         * Ax**2 + Bxy + Cy**2 + Dx + Ey = 1
         *
         * Given a defined (x, y) = (x_n, y_n), each equation resolves to a
         * linear equation (and system homogeneity enables us to normalize the
         * unshown constant onto the right side to correlate with the canonical
         * ellipse equation).
         */
        solVector = solveNormalizedSystem([       
            [points[0][0] * points[0][0], points[0][0] * points[0][1],
                points[0][1] * points[0][1], points[0][0], points[0][1]],
            [points[1][0] * points[1][0], points[1][0] * points[1][1],
                points[1][1] * points[1][1], points[1][0], points[1][1]],
            [points[2][0] * points[2][0], points[2][0] * points[2][1],
                points[2][1] * points[2][1], points[2][0], points[2][1]],
            [points[3][0] * points[3][0], points[3][0] * points[3][1],
                points[3][1] * points[3][1], points[3][0], points[3][1]],
            [points[4][0] * points[4][0], points[4][0] * points[4][1],
                points[4][1] * points[4][1], points[4][0], points[4][1]]
        ], 5);

        if (!solVector) {
            return null;
        }
        
        /* Extract the coefficients from the solution vector. */
        a = solVector[0];  // x ** 2 coefficient
        b = solVector[1];  // x * y coefficient
        c = solVector[2];  // y ** 2 coefficient
        d = solVector[3];  // x coefficient
        e = solVector[4];  // y coefficient
        
        /* Quickly test whether the object is an ellipse or degenerate
         * variant, based on the discriminant. There's a more thorough
         * test for degenerates, too. (See Wikipedia's Ellipse article.)
         * However, I do not deem it optimal here since that would be
         * the minority of applied cases.
         */
        if (!a || !c || b * b - 4.0 * a * c >= 0) {
            return null;
        }

        /* If B != 0, then the ellipse in question is not in its canonical
         * orientation; that is, it is rotated. The canonical equation of
         * the ellipse can be achieved from a rotated coordinated system,
         * corresponding to the equation A'x'**2 + C'y'**2 + D'x' + E'y' = 1.
         * From right-angle trigonometry, the following transforms are
         * derived for all (x, y) in the original coordinate system and for
         * all (x', y') in the rotated coordinate system:
         *
         * x = x' * cos(theta) - y' * sin(theta)
         * y = x' * sin(theta) + y' * cos(theta)
         *
         * When one returns to the original equation with the original five
         * coefficients, x and y may be substituted with the above. In the
         * transformed equation, the x'y' term has an implied zero coefficient.
         * (That is, B' = 0.) When the x'y' terms are combined, the following
         * equation is derived.
         *
         * 0 = -2 * A * sin(theta) * cos(theta) + 2 * C * sin(theta) +
         *     cos(theta) + B * (cos(theta))**2 - B * (sin(theta))**2
         *   = (C - A) * sin(2 * theta) + B * cos(2 * theta)
         * => (A - C) * sin(2 * theta) = B * cos(2 * theta)
         * => tan(2 * theta) = B / (A - C)
         *
         * Therefore:
         * theta = atan(b / (a - c)) / 2 for a - c != 0.
         *
         * If a - c = 0, we take the principal limit of atan(x) / 2
         * as x -> Infinity and find theta = PI / 4.
         */
        theta = ellipseProp.theta =
            (a - c) ? Math.atan(b / (a - c)) / 2.0 : Math.PI / 4.0;

        /* Use right-triangle derived-identities relating the arctangent to
         * the sine and cosine of theta, and combine with the half-angle
         * identities. It probably does just as well just to call
         * Math.cos() and Math.sin(), but this looks cooler (and is
         * theoretically a bit more precise, in particular for the squares).
         */
        hypotenuse = Math.sqrt((a - c) * (a - c) + b * b);
        tmp = Math.abs((a - c) / (2.0 * hypotenuse));  // cos(theta * 2) / 2
        cosThetaSquared = 0.5 + tmp;
        cosTheta = Math.sqrt(cosThetaSquared);
        sinThetaSquared = 0.5 - tmp;
        sinTheta = Math.sqrt(sinThetaSquared);
        if (theta < 0) {
            sinTheta = -sinTheta;
        }

        /* With the aforementioned coordinate system transform, we are now
         * poised to find the major and minor axis lengths for the ellipse in
         * question. We first determine the values of A', D', C', and E' in
         * the above equation. This can be straightforwardly done by combining
         * like terms and comparing coefficients as before.
         */
        tmp = b * cosTheta * sinTheta;  // B * sin(theta) * cos(theta)
        aPrime = a * cosThetaSquared + c * sinThetaSquared + tmp;
        dPrime = d * cosTheta + e * sinTheta;
        cPrime = a * sinThetaSquared + c * cosThetaSquared - tmp;
        ePrime = -d * sinTheta + e * cosTheta;

        /* We are now poised to calculate the ellipse's canonical form. We
         * complete the squares to calculate the offsets and then convert
         * back to the previous coordinate system.
         */
        xOriginPrime = -dPrime / (2.0 * aPrime);
        yOriginPrime = -ePrime / (2.0 * cPrime);
        ellipseProp.xOrigin = xOriginPrime * cosTheta -
                              yOriginPrime * sinTheta;
        ellipseProp.yOrigin = xOriginPrime * sinTheta +
                              yOriginPrime * cosTheta;

        /* We now balance the equation with the squares of the origins and
         * divide. The coefficients of the x ** 2 and y ** 2 terms are the
         * inverse of a ** 2 and b ** 2 in the canonical ellipse formula. The
         * principal square roots are therefore half of the width and half of
         * the height, respectively. xOrigin ** 2 * aPrime and
         * yOrigin ** 2 * cPrime are expanded here for precision.
         */
        normalizer = 1.0 + dPrime * dPrime / (4.0 * aPrime) + 
                           ePrime * ePrime / (4.0 * cPrime);
        ellipseProp.halfWidth = Math.sqrt(normalizer / aPrime);
        ellipseProp.halfHeight = Math.sqrt(normalizer / cPrime);
        ellipseProp.width = ellipseProp.halfWidth * 2.0;
        ellipseProp.height = ellipseProp.halfHeight * 2.0;

        return ellipseProp;
    }

    /**
     * Executes the program and either draws an ellipse passing through five
     * anchor points of a highlighted path, draws a circle passing through
     * three points, or invokes an alert dialog should
     * the conditions not be met or an ellipse not found.
     */
    function main() {
        var currDoc = app.activeDocument,
            inputPath = currDoc.selection[0],
            currLayer = currDoc.activeLayer,
            ellipse,
            ellipseProp,
            circleProp,
            anchors = [],
            i;
            
        if (inputPath) {
            for (i = 0; i < inputPath.pathPoints.length; i += 1) {
                anchors.push(inputPath.pathPoints[i].anchor);
            }
    
            if (anchors.length === 5) {
                ellipseProp = calculateEllipse(anchors);
                
                if (ellipseProp) {
                    // Draw the matching ellipse.
                    ellipse = currLayer.pathItems.ellipse(
                        ellipseProp.yOrigin + ellipseProp.halfHeight,
                        ellipseProp.xOrigin - ellipseProp.halfWidth,
                        ellipseProp.width,
                        ellipseProp.height
                    );
                    
                    // Rotate according to the earlier-calculated rotation
                    // angle.
                    ellipse.rotate(ellipseProp.theta * 180.0 / Math.PI);
                    
                    // Delete original path.
                    inputPath.remove();
                } else {
                    // No ellipse found.
                    alert('No ellipse exists passing through these points.');
                }
            } else if (anchors.length === 3) {
                circleProp = calculateCircle(anchors);
                
                if (circleProp) {
                    // Because circles are ellipses, too!
                    ellipse = currLayer.pathItems.ellipse(
                        circleProp.yOrigin + circleProp.radius,
                        circleProp.xOrigin - circleProp.radius,
                        circleProp.diameter,
                        circleProp.diameter
                    );

                    inputPath.remove();
                } else {
                    alert('No circle exists passing through these points.');
                }
            } else {
                alert('Exactly three or five anchor points must on the'  +
                      'selected path.');
            }
        } else {
            alert('To run, please draw and select a path with exactly three ' +
                  'or five anchor points.');
        }
    }
    
    main();
}());
