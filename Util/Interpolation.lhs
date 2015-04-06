\begin{code}


module Util.Interpolation 
(linear
,linearFit
,lagrange
,Spline
,evalSpline
,linearSpline
,quadraticSpline
,cubicSpline
,flatCubicSpline
,cubicFD
,Poly2
,evalPoly2
,trimPoly2
,biLinear
,biCubic
)where

import Util.Base

import qualified Math.Polynomial as P
import Control.Monad
import Data.List
import Data.Array.Unboxed as V 


\end{code}


1 dimensional interpolation, done with lists of tuples

linear
\begin{code}
linearSlope :: Fractional a => (a,a) -> (a,a) -> a
linearSlope (x,y) (u,v) =  (y - v)/(x - u)

linearExtrapolation :: Num a => (a,a) -> a -> a -> a
linearExtrapolation (x,y) slope x' = y + (x' - x)*slope

linearIntrp :: Fractional a => (a,a) -> (a,a) -> a -> a
linearIntrp pnt1 pnt2 x' = linearExtrapolation pnt1 (linearSlope pnt1 pnt2 ) x'

linear :: (Fractional a,Eq a) => (a,a) -> (a,a) -> P.Poly a
linear (u,v) y = 
	let 
		a = linearSlope (u,v) y 
		b = v - a*u
	in P.addPoly (P.scalePoly a P.x) (P.constPoly b)
	

--
--linear fit by least square
--linearFit :: (Eq a , Fractional a) => [(a,a)] -> P.Poly a
linearFit :: (Eq a , Fractional a) => [(a,a)] -> P.Poly a
linearFit points = 
	let 
		(x',y') = centerFit points
		a = (sum . map (\(u,v) -> (u - x')*(v - y') ) $ points) / (sum . map (\u -> (u - x')^2 ) . map snd $ points )
		b = y' - a*x'
	in P.addPoly (P.scalePoly a P.x) (P.constPoly b)
	

centerFit :: Fractional a => [(a,a)] -> (a,a)
centerFit points = (  (sum (map fst points)) / ( fromIntegral $ length points) , (sum (map snd points)) / (fromIntegral $ length points) )
--
\end{code}



n-th degree polynomials for (n+1) points, lagrange method
NB use only for small n, large n leads to extreme behavior.
\begin{code}

--note if there are any duplicate x coordinates there will be a divide by zero.
lagrange :: (Eq a, Fractional a) => [(a,a)] -> P.Poly a
lagrange points = P.sumPolys . snd $ foldr f (0,[]) points
	where
		--f :: Fractional a => (a,a) -> (Int,[P.Poly a]) -> (Int,[P.Poly a])
		f (x,y) (counter,polys) = ( counter + 1 , (P.scalePoly y $ foldr (g x) P.one (removeNthElement counter points ) ) : polys )
		--g :: Fractional a => a -> (a,a) -> P.Poly a -> P.Poly a
		g x (u,v) p = P.multPoly p $ (P.addPoly P.x $ P.constPoly (-u) ) `P.quotPoly` (P.constPoly (x - u))
--}

--zero indexed
removeNthElement :: Int -> [a] -> [a]
removeNthElement n list = snd $ foldr (\x (m,list') -> if m == n then (m+1,list') else (m+1, (x : list') ) ) (0,[]) list

--one indexed
removeNthElement1 :: Int -> [a] -> [a]
removeNthElement1 n list = snd $ foldr (\x (m,list') -> if m == n then (m+1,list') else (m+1, (x : list') ) ) (1,[]) list

\end{code}

n-th degree polynomial, newton method
\begin{code}

{- This should be equivilent to lagrange, though allowing us to more easily add new points.
newton :: (Eq a, Fractional a) => [(a,a)] -> P.Poly a
newton points =  foldSeq P.addPoly P.zero (\n -> P.scalePoly (a n) $ foldr P.multPoly P.one $ replicate n P.x ) 0 ((length points) - 1)
	where 
		f n = snd $ points !! n
		x n = fst $ points !! n
		a n = (f n - ( sumSeq (\i -> (a i) * prodSeq (\j -> (x n) - (x j) ) 1 (i - 1) ) ) 1 (n - 1) ) / (prodSeq (\j -> (x n) - (x j) ) 1 (n - 1) )
-}

--newtonAddPoint


\end{code}



Polynomial Splines
\begin{code}

--the function returns true if a is in the domain of the associated poly
--i.e., true if we can use the poly to evalute the spline
type Spline a = [ (a -> Bool,P.Poly a) ]


evalSpline :: (MonadPlus m, Num a, Eq a) => Spline a -> a -> m a
evalSpline s x = foldr (\(test,p) y -> if test x then y `mplus` (return $ P.evalPoly p x) else y ) mzero s

splineToPoly :: MonadPlus m => Spline a -> a -> m (P.Poly a)
splineToPoly s x = foldr (\(test,p) y -> if test x then y `mplus` (return p) else y ) mzero s

\end{code}

Note all of these splines operate under the assumption that the points have been ordered with repsect to x 

Linear
Continuous
\begin{code}

containmentTest :: Ord a => (a,a) -> (a,a) -> (a -> Bool)
containmentTest (x,y) (u,v) = \t -> (x <= t) && (t <= u)

--note this currently does not extrapolate beyond the bounds of 
linearSpline :: (Ord a , Fractional a) => [(a,a)] -> Spline a--foldl' due to the importance of ordering.
linearSpline (x:xs) = snd $ foldl' f (x,[]) xs 
	where
		--f :: ((a,a),Spline a) -> (a,a) -> ((a,a),Spline a)
		f  (u,s) v = ( v ,(containmentTest u v , linear u v ) : s)
linearSpline _ = [] 

\end{code}


Naive Quadratic
Continuous first derivatives
\begin{code}

quadraticSpline :: (Ord a , Fractional a) => [(a,a)] -> Spline a
quadraticSpline (x:y:ys) = snd $ foldl' f (y,[(containmentTest x y , linear x y)]) ys
	where
		--f::((a,a),Spline a) -> (a,a) -> ((a,a),Spline a)
		f (u,s) v = ( v ,(containmentTest u v , p ) : s)
			where --see below 
				dsu = ( P.evalPoly (P.polyDeriv . snd . head $ s) (fst u) ) 
				a = ((fst u)*dsu - dsu*(fst v) - (snd u) + (snd v))/((fst u)-(fst v))^2
				b = (-(dsu)*(fst u)^2 +2*(fst u)*(snd u) - 2*(fst u)*(snd v) + dsu*(fst v)^2)/((fst u)-(fst v))^2
				c = (dsu*(fst v)*(fst u)^2+(snd v)*(fst u)^2-(fst u)*dsu*(fst v)^2 - 2*(fst u)*(snd u)*(fst v) + (snd u)*(fst v)^2)/((fst u)-(fst v))^2
				p = P.sumPolys [  P.scalePoly a $ P.multPoly P.x P.x  , P.scalePoly b P.x ,  P.constPoly c ]
quadraticSpline _ = []
				
{-- dsu = derivative of of the previous qudratic when evaluated at u
2*a(fst u) + b = dsu
a(fst u)^2 + b(fst u) + c = (snd u)
a(fst v)^2 + b(fst v) + c = (snd v)
=> By wolfram alpha
Where (fst u) /= (fst v)  
a = ((fst u)*dsu - dsu*(fst v) - (snd u) + (snd v))/((fst u)-(fst v))^2,   
b = (-(dsu)*(fst u)^2 +2*(fst u)*(snd u) - 2*(fst u)*(snd v) + dsu*(fst v)^2)/((fst u)-(fst v))^2,   
c = (dsu*(fst v)*(fst u)^2+(snd v)*(fst u)^2-(fst u)*dsu*(fst v)^2 - 2*(fst u)*(snd u)*(fst v) + (snd u)*(fst v)^2)/((fst u)-(fst v))^2
=> 
--}		




\end{code}


Naive Cubic
Continuous second derivatives after the first point
Note, this tends to make extreme curves for repeated flat points. 
\begin{code}


cubicSpline :: (Ord a , Fractional a) => [(a,a)] -> Spline a
cubicSpline (x:y:z:zs) = snd $ foldl' f (z,quadraticSpline [x,y,z]) zs
	where --f::((a,a),Spline a) -> (a,a) -> ((a,a),Spline a)
		f (u,s) v = ( v ,(containmentTest u v , p ) : s)
			where --see below 
				dsu = ( P.evalPoly (P.polyDeriv . snd . head $ s) (fst u) ) 
				ddsu = ( P.evalPoly (P.polyDeriv . P.polyDeriv . snd . head $ s) (fst u) ) 
				a = ( (-2)*(snd v) + 2*(snd u) + (ddsu)*(fst u)^2 +( (-2)*(ddsu)*(fst v) - 2*(dsu))*(fst u) + (ddsu)*(fst v)^2 + 2*(dsu)*(fst v)) /(2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
				b = -((fst u)*((-6)*(snd v) + 6*(snd u) + 6*(dsu)*(fst v)) + 2*(ddsu)*(fst u)^3 + ((-3)*(ddsu)*(fst v) - 6*(dsu))*(fst u)^2 + (ddsu)*(fst v)^3) /(2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
				c = (((fst u)^2)*((-6)*(snd v) + 6*(snd u) - 3*(ddsu)*(fst v)^2 ) + (ddsu)*(fst u)^4 - 4*(dsu)*(fst u)^3 + (2*(ddsu)*(fst v)^3 + 6*(dsu)*(fst v)^2)*(fst u) - 2*(dsu)*(fst v)^3 ) /( 2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
				d = -(((fst u)^3)*((-2)*(snd v) - 2*(ddsu)*(fst v)^2 - 4*(dsu)*(fst v)) + (fst u)*( (-6)*((fst v)^2)*(snd u) - 2*(dsu)*(fst v)^3 ) + ((fst u)^2)*( 6*(fst v)*(snd u) + (ddsu)*(fst v)^3 + 6*(dsu)*(fst v)^2 ) + 2*((fst v)^3)*(snd u) + (ddsu)*(fst v)*(fst u)^4) /( 2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
				p = P.sumPolys [ P.scalePoly a $ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly b $ P.multPoly P.x P.x  , P.scalePoly c P.x ,  P.constPoly d ] --}
cubicSpline x = quadraticSpline x

{--
As above, dsu = derivative of s at u, ddsu = second derivative of s at u

6a(fst u) + 2b = ddsu
3a(fst u)^2 + 2b(fst u) + c = dsu
a(fst u)^3 + b(fst u)^2 + c(fst u) + d = (snd u)
a(fst v)^3 + b(fst v)^2 + c(fst v) + d = (snd v)

--solving by computer gives,

6*a*(x) + 2*b = p
3*a*(x)^2 + 2*b*(x) + c = q
a*(x)^3 + b*(x)^2 + c*(x) + d = (y) 
a*(v)^3 + b*(v)^2 + c*(v) + d = (z)
=>
a = ( (-2)*z + 2*y + p*x^2 +( (-2)*p*v - 2*q)*x + p*v^2 + 2*q*v) /(2*x^3 - 6*v*x^2 + 6*(v^2)*x - 2*v^3)
b = -(x*((-6)*z + 6*y + 6*q*v) + 2*p*x^3 + ((-3)*p*v - 6*q)*x^2 + p*v^3) /(2*x^3 - 6*v*x^2 + 6*(v^2)*x - 2*v^3)
c = ((x^2)*((-6)*z + 6*y - 3*p*v^2 ) + p*x^4 - 4*q*x^3 + (2*p*v^3 + 6*q*v^2)*x - 2*q*v^3 ) /( 2*x^3 - 6*v*x^2 + 6*(v^2)*x - 2*v^3)
d = -((x^3)*((-2)*z - 2*p*v^2 - 4*q*v) + x*( (-6)*(v^2)*y - 2*q*v^3 ) + (x^2)*( 6*v*y + p*v^3 + 6*q*v^2 ) + 2*(v^3)*y + p*v*x^4) /( 2*x^3 - 6*v*x^2 + 6*(v^2)*x - 2*v^3)
=>
a = ( (-2)*(snd v) + 2*(snd u) + (ddsu)*(fst u)^2 +( (-2)*(ddsu)*(fst v) - 2*(dsu))*(fst u) + (ddsu)*(fst v)^2 + 2*(dsu)*(fst v)) /(2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
b = -((fst u)*((-6)*(snd v) + 6*(snd u) + 6*(dsu)*(fst v)) + 2*(ddsu)*(fst u)^3 + ((-3)*(ddsu)*(fst v) - 6*(dsu))*(fst u)^2 + (ddsu)*(fst v)^3) /(2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
c = (((fst u)^2)*((-6)*(snd v) + 6*(snd u) - 3*(ddsu)*(fst v)^2 ) + (ddsu)*(fst u)^4 - 4*(dsu)*(fst u)^3 + (2*(ddsu)*(fst v)^3 + 6*(dsu)*(fst v)^2)*(fst u) - 2*(dsu)*(fst v)^3 ) /( 2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)
d = -(((fst u)^3)*((-2)*(snd v) - 2*(ddsu)*(fst v)^2 - 4*(dsu)*(fst v)) + (fst u)*( (-6)*((fst v)^2)*(snd u) - 2*(dsu)*(fst v)^3 ) + ((fst u)^2)*( 6*(fst v)*(snd u) + (ddsu)*(fst v)^3 + 6*(dsu)*(fst v)^2 ) + 2*((fst v)^3)*(snd u) + (ddsu)*(fst v)*(fst u)^4) /( 2*(fst u)^3 - 6*(fst v)*(fst u)^2 + 6*((fst v)^2)*(fst u) - 2*(fst v)^3)

--}


\end{code}

FlatCubic
Here we assume a zero derivative at every given point
\begin{code}

flatCubicSpline :: (Ord a , Fractional a) => [(a,a)] -> Spline a
flatCubicSpline (x:y:ys) = snd $ foldl' f (x,[]) (y:ys)
	where --f::((a,a),Spline a) -> (a,a) -> ((a,a),Spline a)
		f (u,s) v = ( v ,(containmentTest u v , p ) : s)
			where --from above
				a =  -(2*(snd u)-2*(snd v))/(-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3)
				b = ((3*(snd u) - 3*(snd v))*(fst v) + (3*(snd u) - 3*(snd v))*(fst u)) / (-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3) 
				c = -(6*(snd u) - 6*(snd v))*(fst u)*(fst v) /(-(fst v)^3 + 3*(fst u)*((fst v)^2) - 3*((fst u)^2)*(fst v) + (fst u)^3)
				d = ((-(snd u))*(fst v)^3 + 3*(snd u)*(fst u)*(fst v)^2 - 3*(snd v)*((fst u)^2)*(fst v) + (snd v)*(fst u)^3 )/(-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3)
				p = P.sumPolys [ P.scalePoly a $ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly b $ P.multPoly P.x P.x  , P.scalePoly c P.x ,  P.constPoly d ] --}
flatCubicSpline x = quadraticSpline x

{-

3a(fst u)^2 + 2b(fst u) + c = 0
3a(fst v)^2 + 2b(fst v) + c = 0
a(fst u)^3 + b(fst u)^2 + c(fst u) + d = (snd u)
a(fst v)^3 + b(fst v)^2 + c(fst v) + d = (snd v)

3a(u)^2 + 2b(u) + c = 0
3a(v)^2 + 2b(v) + c = 0
a(u)^3 + b(u)^2 + c(u) + d = (t)
a(v)^3 + b(v)^2 + c(v) + d = (r)

a =  -(2*(snd u)-2*(snd v))/(-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3), 
b = ((3*(snd u) - 3*(snd v))*(fst v) + (3*(snd u) - 3*(snd v))*(fst u)) / (-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3), 
c = -(6*(snd u) - 6*(snd v))*(fst u)*(fst v) /(-(fst v)^3 + 3*(fst u)*((fst v)^2) - 3*((fst u)^2)*(fst v) + (fst u)^3), 
d = ((-(snd u))*(fst v)^3 + 3*(snd u)*(fst u)*(fst v)^2 - 3*(snd v)*((fst u)^2)*(fst v) + (snd v)*(fst u)^3 )/(-(fst v)^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3)


-}


\end{code}




Cubic - finite differance tangents
Note that because of the uniqueness of solution, this should be the same as hermite with finite differances.
\begin{code}

cubicFD :: (Ord a , Fractional a) => [(a,a)] -> Spline a
cubicFD points@(t:r:s:ys) = (\(_,_,_,s') -> s' ) $ foldl' f (head points',t,r,[]) ( tail . tail . tail $ points')
	where
		points' = appendPoints points
		f (x,u,v,s) y = ( u, v , y ,(containmentTest u v , p ) : s)
			where
				dsu = ( (snd u) - (snd x) ) / (2*( (fst u) - (fst x) )) + ( (snd v) - (snd u) ) / (2*( (fst v) -(fst u) )) 
				dsv = ( (snd y) - (snd v) ) / (2*( (fst y) - (fst v) )) + ( (snd v) - (snd u) ) / (2*( (fst v) -(fst u) )) 
				a = ((-(fst v))*(dsv) + (fst u)*(dsv) + ( (fst u) - (fst v) )*(dsu) + 2*(snd v) - 2*(snd u) ) / ( (-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )
				b = -( (fst u)*( (-(fst v))*(dsv) + 3*(snd v) - 3*(snd u) ) - ((fst v)^2)*(dsv) + 2*((fst u)^2)*(dsv) + ((-2)*(fst v)^2 + (fst u)*(fst v) + (fst u)^2 )*(dsu) + (3*(snd v) - 3*(snd u))*(fst v) ) / ((-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )
				c = ((fst u)*( ( 6*(snd v) - 6*(snd u) )*(fst v) - 2*((fst v)^2)*(dsv) ) + ((fst u)^2)*(fst v)*(dsv) + ((fst u)^3)*(dsv) + ((-(fst v))^3 - (fst u)*(fst v)^2 + 2*((fst u)^2)*(fst v) )*(dsu) ) / ((-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3)
				d = -(((fst u)^2)*(3*(snd v)*(fst v) - ((fst v)^2)*(dsv)) + ((fst u)^3)*( (fst v)*(dsv) - (snd v) ) + ( ((fst u)^2)*((fst v)^2) - (fst u)*(fst v)^3 )*(dsu) + (snd u)*(fst v)^3 - 3*(snd u)*(fst u)*(fst v)^2 ) / ( (-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )
				p = P.sumPolys [ P.scalePoly a $ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly b $ P.multPoly P.x P.x  , P.scalePoly c P.x ,  P.constPoly d ] --}
--

appendPoints :: Fractional a => [(a,a)] -> [(a,a)]
appendPoints xs@(x:y:ys) = 
	let 
		u = last xs
		v = last $ init xs
	in [ ( (fst x) - 1 , linearIntrp x y ((fst x) - 1) ) ] 
		++ xs 
		++ [( (fst u) + 1 , linearIntrp u v ((fst u) + 1)) ] 
appendPoints _ = [] 


{-


3a(fst u)^2 + 2b(fst u) + c = dsu
3a(fst v)^2 + 2b(fst v) + c = dsv
a(fst u)^3 + b(fst u)^2 + c(fst u) + d = (snd u)
a(fst v)^3 + b(fst v)^2 + c(fst v) + d = (snd v) 

3*a*(u)^2 + 2*b*(u) + c = w
3*a*(v)^2 + 2*b*(v) + c = z
a*(u)^3 + b*(u)^2 + c*(u) + d = q
a*(v)^3 + b*(v)^2 + c*(v) + d = r

a = ((-v)*z + u*z + ( u - v )*w + 2*r - 2*q ) / ( (-v)^3 + 3*u*v^2 - 3*(u^2)*v + u^3 )
b = -( u*( (-v)*z + 3*r - 3*q ) - (v^2)*z + 2*(u^2)*z + ((-2)*v^2 + u*v + u^2 )*w + (3*r - 3*q)*v ) / ((-v)^3 + 3*u*v^2 - 3*(u^2)*v + u^3 )
c = (u*( ( 6*r - 6*q )*v - 2*(v^2)*z ) + (u^2)*v*z + (u^3)*z + ((-v)^3 - u*v^2 + 2*(u^2)*v )*w ) / ((-v)^3 + 3*u*v^2 - 3*(u^2)*v + u^3)
d = -((u^2)*(3*r*v - (v^2)*z) + (u^3)*( v*z - r ) + ( (u^2)*(v^2) - u*v^3 )*w + q*v^3 - 3*q*u*v^2 ) / ( (-v)^3 + 3*u*v^2 - 3*(u^2)*v + u^3 )

a = ((-(fst v))*(dsv) + (fst u)*(dsv) + ( (fst u) - (fst v) )*(dsu) + 2*(snd v) - 2*(snd u) ) / ( (-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )
b = -( (fst u)*( (-(fst v))*(dsv) + 3*(snd v) - 3*(snd u) ) - ((fst v)^2)*(dsv) + 2*((fst u)^2)*(dsv) + ((-2)*(fst v)^2 + (fst u)*(fst v) + (fst u)^2 )*(dsu) + (3*(snd v) - 3*(snd u))*(fst v) ) / ((-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )
c = ((fst u)*( ( 6*(snd v) - 6*(snd u) )*(fst v) - 2*((fst v)^2)*(dsv) ) + ((fst u)^2)*(fst v)*(dsv) + ((fst u)^3)*(dsv) + ((-(fst v))^3 - (fst u)*(fst v)^2 + 2*((fst u)^2)*(fst v) )*(dsu) ) / ((-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3)
d = -(((fst u)^2)*(3*(snd v)*(fst v) - ((fst v)^2)*(dsv)) + ((fst u)^3)*( (fst v)*(dsv) - (snd v) ) + ( ((fst u)^2)*((fst v)^2) - (fst u)*(fst v)^3 )*(dsu) + (snd u)*(fst v)^3 - 3*(snd u)*(fst u)*(fst v)^2 ) / ( (-(fst v))^3 + 3*(fst u)*(fst v)^2 - 3*((fst u)^2)*(fst v) + (fst u)^3 )

-}

\end{code}








2 dimensional interpolation:

Rough 2d polys
\begin{code}

-- [x^0[y^0,y^1 ...] , x^1[y^0,y^1 ...]  ... ]
newtype Poly2 a = Poly2 [[a]] deriving (Eq,Show,Ord)

evalPoly2 :: Num a => Poly2 a -> (a,a) -> a
evalPoly2 (Poly2 list) (x,y) = foldr f 0 [0..n]
	where
		n = (length list) - 1
		f i total' = foldr g total' [0..m]
			where
				m = length (list !! i) - 1
				g j total = total + ((list !! i) !! j)*(x^i)*(y^j)


trimPoly2 :: (Num a,Eq a) => Poly2 a -> Poly2 a
trimPoly2 (Poly2 list) = Poly2 $ map (reverse . dropWhile (0 ==) . reverse ) list



--3 + y^2 -y^3 + xy + 4xy^3 - x^3
t1 = Poly2 [[3,0,1,-1],[0,1,0,4],[0,0,0,0],[-1,0,0,0]]
{-t2
t3
t4
-}

\end{code}





--note we assume that arrays are indexed from the bottom left with (x,y)
i.e.
a_01 a_11
a_00 a_10
--We also assume that every point in the array has value, with the distance between points being scaled to 1.

Bilinear
\begin{code}
biLinear :: RealFrac a => Array (Int,Int) a -> (a,a) -> Poly2 a
biLinear arr (x,y) = Poly2 [[a00,a01],[a10,a11]]
	where
		getPoint (u,v) = if and [ n_ <= u , u <= n' , m_ <= v , v <= m' ] then arr V.! (u,v) else 0
			where 
				((n_,m_),(n',m')) = V.bounds arr
		x0 = fromIntegral $ floor x
		x1 = fromIntegral $ 1 + floor x
		y0 = fromIntegral $ floor y
		y1 = fromIntegral $ 1 + floor y
		p00 = getPoint (floor x,floor y)
		p01 = getPoint (floor x,1 + floor y)
		p10 = getPoint (1 + floor x,floor y)
		p11 = getPoint (1 + floor x,1 + floor y)
		a00 = (p00*x1*y1-p01*x1*y0-p10*x0*y1+p11*x0*y0)/((x0-x1)*(y0-y1))
		a01 = -(p00*x1-p01*x1-p10*x0+p11*x0)/((x0-x1)*(y0-y1))
		a10 = -(p00*y1-p01*y0-p10*y1+p11*y0)/((x0-x1)*(y0-y1))
		a11 = (p00-p01-p10+p11)/(x0*y0-x0*y1-x1*y0+x1*y1)
--

{--
Of form
a00 + a01y + a10x + a11xy

Where we want to garruntee
{a00 + a01*(y0) + a10*(x0) + a11*(x0)*(y0) = p00,
a00 + a01*(y1) + a10*(x0) + a11*(x0)*(y1) = p01,
a00 + a01*(y0) + a10*(x1) + a11*(x1)*(y0) = p10,
a00 + a01*(y1) + a10*(x1) + a11*(x1)*(y1) = p11}

Solving by maple for {a00,a01,a10,a11} gives
a00 = (p00*x1*y1-p01*x1*y0-p10*x0*y1+p11*x0*y0)/((x0-x1)*(y0-y1))
a01 = -(p00*x1-p01*x1-p10*x0+p11*x0)/((x0-x1)*(y0-y1))
a10 = -(p00*y1-p01*y0-p10*y1+p11*y0)/((x0-x1)*(y0-y1))
a11 = (p00-p01-p10+p11)/(x0*y0-x0*y1-x1*y0+x1*y1)

--}

\end{code}











Bicubic interpolation; done with Arrays
Based on http://en.wikipedia.org/wiki/Bicubic_interpolation

\begin{code}
--We assume 1 between each point
biCubic :: RealFrac a => Array (Int,Int) a -> (a,a) -> Poly2 a
biCubic arr (x,y) = Poly2 [[a00,a01,a02,a03],[a10,a11,a12,a13],[a20,a21,a22,a23],[a30,a31,a32,a33]]
	where
		getPoint (u,v) = if and [ n_ <= u , u <= n' , m_ <= v , v <= m' ] then arr V.! (u,v) else 0
			where 
				((n_,m_),(n',m')) = V.bounds arr
		--The 4 surronding coordinates
		x0 = fromIntegral $ floor x
		x1 = fromIntegral $ 1 + floor x
		y0 = fromIntegral $ floor y
		y1 = fromIntegral $ 1 + floor y
		--The 16 Points we will be using
		p00 = getPoint (floor x,floor y)
		p01 = getPoint (floor x,1 + floor y)
		p10 = getPoint (1 + floor x,floor y)
		p11 = getPoint (1 + floor x,1 + floor y)
		p02 = getPoint (floor x,2 + floor y)
		p0_1 = getPoint (floor x, -1 + floor y)
		p20 =  getPoint (2 + floor x,floor y)
		p_10 = getPoint ( -1 + floor x,floor y)
		p12 = getPoint (1 + floor x,2 + floor y)
		p1_1 = getPoint (1+ floor x, -1 + floor y)
		p21 =  getPoint (2 + floor x,1 + floor y)
		p_11 = getPoint ( -1 + floor x,1 + floor y)
		p22 = getPoint (2 + floor x,2 + floor y)
		p2_1 = getPoint (2 + floor x,-1 + floor y)
		p_12 = getPoint (-1 + floor x,2 + floor y)
		p_1_1 = getPoint (-1 + floor x,-1 + floor y)
		--finite differance derivatives
		dx00 = (p10 - p_10)/2
		dx01 = (p11 - p_11)/2
		dx10 = (p20 - p00)/2
		dx11 = (p21 - p01)/2
		dy00 = (p01 - p0_1)/2
		dy01 = (p02 - p00)/2
		dy10 = (p11 - p1_1)/2
		dy11 = (p12 - p10)/2
		dxy00 = (p11 - p1_1 - p_11 + p_1_1)/4
		dxy01 = (p12 - p10 - p_12 + p_10)/4
		dxy10 = (p21 - p2_1 - p01 + p0_1)/4
		dxy11 = (p22 - p20 - p02 + p00)/4
		--from Maple
		a00 = -(-dxy00*x0^2*x1^2*y0^2*y1^2+dxy00*x0^2*x1^2*y0*y1^3+dxy00*x0*x1^3*y0^2*y1^2-dxy00*x0*x1^3*y0*y1^3-dxy01*x0^2*x1^2*y0^3*y1+dxy01*x0^2*x1^2*y0^2*y1^2+dxy01*x0*x1^3*y0^3*y1-dxy01*x0*x1^3*y0^2*y1^2-dxy10*x0^3*x1*y0^2*y1^2+dxy10*x0^3*x1*y0*y1^3+dxy10*x0^2*x1^2*y0^2*y1^2-dxy10*x0^2*x1^2*y0*y1^3-dxy11*x0^3*x1*y0^3*y1+dxy11*x0^3*x1*y0^2*y1^2+dxy11*x0^2*x1^2*y0^3*y1-dxy11*x0^2*x1^2*y0^2*y1^2+3*dx00*x0^2*x1^2*y0*y1^2-dx00*x0^2*x1^2*y1^3-3*dx00*x0*x1^3*y0*y1^2+dx00*x0*x1^3*y1^3+dx01*x0^2*x1^2*y0^3-3*dx01*x0^2*x1^2*y0^2*y1-dx01*x0*x1^3*y0^3+3*dx01*x0*x1^3*y0^2*y1+3*dx10*x0^3*x1*y0*y1^2-dx10*x0^3*x1*y1^3-3*dx10*x0^2*x1^2*y0*y1^2+dx10*x0^2*x1^2*y1^3+dx11*x0^3*x1*y0^3-3*dx11*x0^3*x1*y0^2*y1-dx11*x0^2*x1^2*y0^3+3*dx11*x0^2*x1^2*y0^2*y1+3*dy00*x0*x1^2*y0^2*y1^2-3*dy00*x0*x1^2*y0*y1^3-dy00*x1^3*y0^2*y1^2+dy00*x1^3*y0*y1^3+3*dy01*x0*x1^2*y0^3*y1-3*dy01*x0*x1^2*y0^2*y1^2-dy01*x1^3*y0^3*y1+dy01*x1^3*y0^2*y1^2+dy10*x0^3*y0^2*y1^2-dy10*x0^3*y0*y1^3-3*dy10*x0^2*x1*y0^2*y1^2+3*dy10*x0^2*x1*y0*y1^3+dy11*x0^3*y0^3*y1-dy11*x0^3*y0^2*y1^2-3*dy11*x0^2*x1*y0^3*y1+3*dy11*x0^2*x1*y0^2*y1^2-9*p00*x0*x1^2*y0*y1^2+3*p00*x0*x1^2*y1^3+3*p00*x1^3*y0*y1^2-p00*x1^3*y1^3-3*p01*x0*x1^2*y0^3+9*p01*x0*x1^2*y0^2*y1+p01*x1^3*y0^3-3*p01*x1^3*y0^2*y1-3*p10*x0^3*y0*y1^2+p10*x0^3*y1^3+9*p10*x0^2*x1*y0*y1^2-3*p10*x0^2*x1*y1^3-p11*x0^3*y0^3+3*p11*x0^3*y0^2*y1+3*p11*x0^2*x1*y0^3-9*p11*x0^2*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
		a01 = (-2*dxy00*x0^2*x1^2*y0^2*y1+dxy00*x0^2*x1^2*y0*y1^2+dxy00*x0^2*x1^2*y1^3+2*dxy00*x0*x1^3*y0^2*y1-dxy00*x0*x1^3*y0*y1^2-dxy00*x0*x1^3*y1^3-dxy01*x0^2*x1^2*y0^3-dxy01*x0^2*x1^2*y0^2*y1+2*dxy01*x0^2*x1^2*y0*y1^2+dxy01*x0*x1^3*y0^3+dxy01*x0*x1^3*y0^2*y1-2*dxy01*x0*x1^3*y0*y1^2-2*dxy10*x0^3*x1*y0^2*y1+dxy10*x0^3*x1*y0*y1^2+dxy10*x0^3*x1*y1^3+2*dxy10*x0^2*x1^2*y0^2*y1-dxy10*x0^2*x1^2*y0*y1^2-dxy10*x0^2*x1^2*y1^3-dxy11*x0^3*x1*y0^3-dxy11*x0^3*x1*y0^2*y1+2*dxy11*x0^3*x1*y0*y1^2+dxy11*x0^2*x1^2*y0^3+dxy11*x0^2*x1^2*y0^2*y1-2*dxy11*x0^2*x1^2*y0*y1^2+6*dx00*x0^2*x1^2*y0*y1-6*dx00*x0*x1^3*y0*y1-6*dx01*x0^2*x1^2*y0*y1+6*dx01*x0*x1^3*y0*y1+6*dx10*x0^3*x1*y0*y1-6*dx10*x0^2*x1^2*y0*y1-6*dx11*x0^3*x1*y0*y1+6*dx11*x0^2*x1^2*y0*y1+6*dy00*x0*x1^2*y0^2*y1-3*dy00*x0*x1^2*y0*y1^2-3*dy00*x0*x1^2*y1^3-2*dy00*x1^3*y0^2*y1+dy00*x1^3*y0*y1^2+dy00*x1^3*y1^3+3*dy01*x0*x1^2*y0^3+3*dy01*x0*x1^2*y0^2*y1-6*dy01*x0*x1^2*y0*y1^2-dy01*x1^3*y0^3-dy01*x1^3*y0^2*y1+2*dy01*x1^3*y0*y1^2+2*dy10*x0^3*y0^2*y1-dy10*x0^3*y0*y1^2-dy10*x0^3*y1^3-6*dy10*x0^2*x1*y0^2*y1+3*dy10*x0^2*x1*y0*y1^2+3*dy10*x0^2*x1*y1^3+dy11*x0^3*y0^3+dy11*x0^3*y0^2*y1-2*dy11*x0^3*y0*y1^2-3*dy11*x0^2*x1*y0^3-3*dy11*x0^2*x1*y0^2*y1+6*dy11*x0^2*x1*y0*y1^2-18*p00*x0*x1^2*y0*y1+6*p00*x1^3*y0*y1+18*p01*x0*x1^2*y0*y1-6*p01*x1^3*y0*y1-6*p10*x0^3*y0*y1+18*p10*x0^2*x1*y0*y1+6*p11*x0^3*y0*y1-18*p11*x0^2*x1*y0*y1)/((x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3)*(x0-x1))
		a02 = -(-dxy00*x0^2*x1^2*y0^2-dxy00*x0^2*x1^2*y0*y1+2*dxy00*x0^2*x1^2*y1^2+dxy00*x0*x1^3*y0^2+dxy00*x0*x1^3*y0*y1-2*dxy00*x0*x1^3*y1^2-2*dxy01*x0^2*x1^2*y0^2+dxy01*x0^2*x1^2*y0*y1+dxy01*x0^2*x1^2*y1^2+2*dxy01*x0*x1^3*y0^2-dxy01*x0*x1^3*y0*y1-dxy01*x0*x1^3*y1^2-dxy10*x0^3*x1*y0^2-dxy10*x0^3*x1*y0*y1+2*dxy10*x0^3*x1*y1^2+dxy10*x0^2*x1^2*y0^2+dxy10*x0^2*x1^2*y0*y1-2*dxy10*x0^2*x1^2*y1^2-2*dxy11*x0^3*x1*y0^2+dxy11*x0^3*x1*y0*y1+dxy11*x0^3*x1*y1^2+2*dxy11*x0^2*x1^2*y0^2-dxy11*x0^2*x1^2*y0*y1-dxy11*x0^2*x1^2*y1^2+3*dx00*x0^2*x1^2*y0+3*dx00*x0^2*x1^2*y1-3*dx00*x0*x1^3*y0-3*dx00*x0*x1^3*y1-3*dx01*x0^2*x1^2*y0-3*dx01*x0^2*x1^2*y1+3*dx01*x0*x1^3*y0+3*dx01*x0*x1^3*y1+3*dx10*x0^3*x1*y0+3*dx10*x0^3*x1*y1-3*dx10*x0^2*x1^2*y0-3*dx10*x0^2*x1^2*y1-3*dx11*x0^3*x1*y0-3*dx11*x0^3*x1*y1+3*dx11*x0^2*x1^2*y0+3*dx11*x0^2*x1^2*y1+3*dy00*x0*x1^2*y0^2+3*dy00*x0*x1^2*y0*y1-6*dy00*x0*x1^2*y1^2-dy00*x1^3*y0^2-dy00*x1^3*y0*y1+2*dy00*x1^3*y1^2+6*dy01*x0*x1^2*y0^2-3*dy01*x0*x1^2*y0*y1-3*dy01*x0*x1^2*y1^2-2*dy01*x1^3*y0^2+dy01*x1^3*y0*y1+dy01*x1^3*y1^2+dy10*x0^3*y0^2+dy10*x0^3*y0*y1-2*dy10*x0^3*y1^2-3*dy10*x0^2*x1*y0^2-3*dy10*x0^2*x1*y0*y1+6*dy10*x0^2*x1*y1^2+2*dy11*x0^3*y0^2-dy11*x0^3*y0*y1-dy11*x0^3*y1^2-6*dy11*x0^2*x1*y0^2+3*dy11*x0^2*x1*y0*y1+3*dy11*x0^2*x1*y1^2-9*p00*x0*x1^2*y0-9*p00*x0*x1^2*y1+3*p00*x1^3*y0+3*p00*x1^3*y1+9*p01*x0*x1^2*y0+9*p01*x0*x1^2*y1-3*p01*x1^3*y0-3*p01*x1^3*y1-3*p10*x0^3*y0-3*p10*x0^3*y1+9*p10*x0^2*x1*y0+9*p10*x0^2*x1*y1+3*p11*x0^3*y0+3*p11*x0^3*y1-9*p11*x0^2*x1*y0-9*p11*x0^2*x1*y1)/((y0-y1)*(x0^3*y0^2-2*x0^3*y0*y1+x0^3*y1^2-3*x0^2*x1*y0^2+6*x0^2*x1*y0*y1-3*x0^2*x1*y1^2+3*x0*x1^2*y0^2-6*x0*x1^2*y0*y1+3*x0*x1^2*y1^2-x1^3*y0^2+2*x1^3*y0*y1-x1^3*y1^2))
		a03 = (-dxy00*x0^2*x1^2*y0+dxy00*x0^2*x1^2*y1+dxy00*x0*x1^3*y0-dxy00*x0*x1^3*y1-dxy01*x0^2*x1^2*y0+dxy01*x0^2*x1^2*y1+dxy01*x0*x1^3*y0-dxy01*x0*x1^3*y1-dxy10*x0^3*x1*y0+dxy10*x0^3*x1*y1+dxy10*x0^2*x1^2*y0-dxy10*x0^2*x1^2*y1-dxy11*x0^3*x1*y0+dxy11*x0^3*x1*y1+dxy11*x0^2*x1^2*y0-dxy11*x0^2*x1^2*y1+2*dx00*x0^2*x1^2-2*dx00*x0*x1^3-2*dx01*x0^2*x1^2+2*dx01*x0*x1^3+2*dx10*x0^3*x1-2*dx10*x0^2*x1^2-2*dx11*x0^3*x1+2*dx11*x0^2*x1^2+3*dy00*x0*x1^2*y0-3*dy00*x0*x1^2*y1-dy00*x1^3*y0+dy00*x1^3*y1+3*dy01*x0*x1^2*y0-3*dy01*x0*x1^2*y1-dy01*x1^3*y0+dy01*x1^3*y1+dy10*x0^3*y0-dy10*x0^3*y1-3*dy10*x0^2*x1*y0+3*dy10*x0^2*x1*y1+dy11*x0^3*y0-dy11*x0^3*y1-3*dy11*x0^2*x1*y0+3*dy11*x0^2*x1*y1-6*p00*x0*x1^2+2*p00*x1^3+6*p01*x0*x1^2-2*p01*x1^3-2*p10*x0^3+6*p10*x0^2*x1+2*p11*x0^3-6*p11*x0^2*x1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
		a10 = (-2*dxy00*x0^2*x1*y0^2*y1^2+2*dxy00*x0^2*x1*y0*y1^3+dxy00*x0*x1^2*y0^2*y1^2-dxy00*x0*x1^2*y0*y1^3+dxy00*x1^3*y0^2*y1^2-dxy00*x1^3*y0*y1^3-2*dxy01*x0^2*x1*y0^3*y1+2*dxy01*x0^2*x1*y0^2*y1^2+dxy01*x0*x1^2*y0^3*y1-dxy01*x0*x1^2*y0^2*y1^2+dxy01*x1^3*y0^3*y1-dxy01*x1^3*y0^2*y1^2-dxy10*x0^3*y0^2*y1^2+dxy10*x0^3*y0*y1^3-dxy10*x0^2*x1*y0^2*y1^2+dxy10*x0^2*x1*y0*y1^3+2*dxy10*x0*x1^2*y0^2*y1^2-2*dxy10*x0*x1^2*y0*y1^3-dxy11*x0^3*y0^3*y1+dxy11*x0^3*y0^2*y1^2-dxy11*x0^2*x1*y0^3*y1+dxy11*x0^2*x1*y0^2*y1^2+2*dxy11*x0*x1^2*y0^3*y1-2*dxy11*x0*x1^2*y0^2*y1^2+6*dx00*x0^2*x1*y0*y1^2-2*dx00*x0^2*x1*y1^3-3*dx00*x0*x1^2*y0*y1^2+dx00*x0*x1^2*y1^3-3*dx00*x1^3*y0*y1^2+dx00*x1^3*y1^3+2*dx01*x0^2*x1*y0^3-6*dx01*x0^2*x1*y0^2*y1-dx01*x0*x1^2*y0^3+3*dx01*x0*x1^2*y0^2*y1-dx01*x1^3*y0^3+3*dx01*x1^3*y0^2*y1+3*dx10*x0^3*y0*y1^2-dx10*x0^3*y1^3+3*dx10*x0^2*x1*y0*y1^2-dx10*x0^2*x1*y1^3-6*dx10*x0*x1^2*y0*y1^2+2*dx10*x0*x1^2*y1^3+dx11*x0^3*y0^3-3*dx11*x0^3*y0^2*y1+dx11*x0^2*x1*y0^3-3*dx11*x0^2*x1*y0^2*y1-2*dx11*x0*x1^2*y0^3+6*dx11*x0*x1^2*y0^2*y1+6*dy00*x0*x1*y0^2*y1^2-6*dy00*x0*x1*y0*y1^3+6*dy01*x0*x1*y0^3*y1-6*dy01*x0*x1*y0^2*y1^2-6*dy10*x0*x1*y0^2*y1^2+6*dy10*x0*x1*y0*y1^3-6*dy11*x0*x1*y0^3*y1+6*dy11*x0*x1*y0^2*y1^2-18*p00*x0*x1*y0*y1^2+6*p00*x0*x1*y1^3-6*p01*x0*x1*y0^3+18*p01*x0*x1*y0^2*y1+18*p10*x0*x1*y0*y1^2-6*p10*x0*x1*y1^3+6*p11*x0*x1*y0^3-18*p11*x0*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
		a11 = -(-4*dxy00*x0^2*x1*y0^2*y1+2*dxy00*x0^2*x1*y0*y1^2+2*dxy00*x0^2*x1*y1^3+2*dxy00*x0*x1^2*y0^2*y1-dxy00*x0*x1^2*y0*y1^2-dxy00*x0*x1^2*y1^3+2*dxy00*x1^3*y0^2*y1-dxy00*x1^3*y0*y1^2-dxy00*x1^3*y1^3-2*dxy01*x0^2*x1*y0^3-2*dxy01*x0^2*x1*y0^2*y1+4*dxy01*x0^2*x1*y0*y1^2+dxy01*x0*x1^2*y0^3+dxy01*x0*x1^2*y0^2*y1-2*dxy01*x0*x1^2*y0*y1^2+dxy01*x1^3*y0^3+dxy01*x1^3*y0^2*y1-2*dxy01*x1^3*y0*y1^2-2*dxy10*x0^3*y0^2*y1+dxy10*x0^3*y0*y1^2+dxy10*x0^3*y1^3-2*dxy10*x0^2*x1*y0^2*y1+dxy10*x0^2*x1*y0*y1^2+dxy10*x0^2*x1*y1^3+4*dxy10*x0*x1^2*y0^2*y1-2*dxy10*x0*x1^2*y0*y1^2-2*dxy10*x0*x1^2*y1^3-dxy11*x0^3*y0^3-dxy11*x0^3*y0^2*y1+2*dxy11*x0^3*y0*y1^2-dxy11*x0^2*x1*y0^3-dxy11*x0^2*x1*y0^2*y1+2*dxy11*x0^2*x1*y0*y1^2+2*dxy11*x0*x1^2*y0^3+2*dxy11*x0*x1^2*y0^2*y1-4*dxy11*x0*x1^2*y0*y1^2+12*dx00*x0^2*x1*y0*y1-6*dx00*x0*x1^2*y0*y1-6*dx00*x1^3*y0*y1-12*dx01*x0^2*x1*y0*y1+6*dx01*x0*x1^2*y0*y1+6*dx01*x1^3*y0*y1+6*dx10*x0^3*y0*y1+6*dx10*x0^2*x1*y0*y1-12*dx10*x0*x1^2*y0*y1-6*dx11*x0^3*y0*y1-6*dx11*x0^2*x1*y0*y1+12*dx11*x0*x1^2*y0*y1+12*dy00*x0*x1*y0^2*y1-6*dy00*x0*x1*y0*y1^2-6*dy00*x0*x1*y1^3+6*dy01*x0*x1*y0^3+6*dy01*x0*x1*y0^2*y1-12*dy01*x0*x1*y0*y1^2-12*dy10*x0*x1*y0^2*y1+6*dy10*x0*x1*y0*y1^2+6*dy10*x0*x1*y1^3-6*dy11*x0*x1*y0^3-6*dy11*x0*x1*y0^2*y1+12*dy11*x0*x1*y0*y1^2-36*p00*x0*x1*y0*y1+36*p01*x0*x1*y0*y1+36*p10*x0*x1*y0*y1-36*p11*x0*x1*y0*y1)/((x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3)*(x0-x1))
		a12 = (-2*dxy00*x0^2*x1*y0^2-2*dxy00*x0^2*x1*y0*y1+4*dxy00*x0^2*x1*y1^2+dxy00*x0*x1^2*y0^2+dxy00*x0*x1^2*y0*y1-2*dxy00*x0*x1^2*y1^2+dxy00*x1^3*y0^2+dxy00*x1^3*y0*y1-2*dxy00*x1^3*y1^2-4*dxy01*x0^2*x1*y0^2+2*dxy01*x0^2*x1*y0*y1+2*dxy01*x0^2*x1*y1^2+2*dxy01*x0*x1^2*y0^2-dxy01*x0*x1^2*y0*y1-dxy01*x0*x1^2*y1^2+2*dxy01*x1^3*y0^2-dxy01*x1^3*y0*y1-dxy01*x1^3*y1^2-dxy10*x0^3*y0^2-dxy10*x0^3*y0*y1+2*dxy10*x0^3*y1^2-dxy10*x0^2*x1*y0^2-dxy10*x0^2*x1*y0*y1+2*dxy10*x0^2*x1*y1^2+2*dxy10*x0*x1^2*y0^2+2*dxy10*x0*x1^2*y0*y1-4*dxy10*x0*x1^2*y1^2-2*dxy11*x0^3*y0^2+dxy11*x0^3*y0*y1+dxy11*x0^3*y1^2-2*dxy11*x0^2*x1*y0^2+dxy11*x0^2*x1*y0*y1+dxy11*x0^2*x1*y1^2+4*dxy11*x0*x1^2*y0^2-2*dxy11*x0*x1^2*y0*y1-2*dxy11*x0*x1^2*y1^2+6*dx00*x0^2*x1*y0+6*dx00*x0^2*x1*y1-3*dx00*x0*x1^2*y0-3*dx00*x0*x1^2*y1-3*dx00*x1^3*y0-3*dx00*x1^3*y1-6*dx01*x0^2*x1*y0-6*dx01*x0^2*x1*y1+3*dx01*x0*x1^2*y0+3*dx01*x0*x1^2*y1+3*dx01*x1^3*y0+3*dx01*x1^3*y1+3*dx10*x0^3*y0+3*dx10*x0^3*y1+3*dx10*x0^2*x1*y0+3*dx10*x0^2*x1*y1-6*dx10*x0*x1^2*y0-6*dx10*x0*x1^2*y1-3*dx11*x0^3*y0-3*dx11*x0^3*y1-3*dx11*x0^2*x1*y0-3*dx11*x0^2*x1*y1+6*dx11*x0*x1^2*y0+6*dx11*x0*x1^2*y1+6*dy00*x0*x1*y0^2+6*dy00*x0*x1*y0*y1-12*dy00*x0*x1*y1^2+12*dy01*x0*x1*y0^2-6*dy01*x0*x1*y0*y1-6*dy01*x0*x1*y1^2-6*dy10*x0*x1*y0^2-6*dy10*x0*x1*y0*y1+12*dy10*x0*x1*y1^2-12*dy11*x0*x1*y0^2+6*dy11*x0*x1*y0*y1+6*dy11*x0*x1*y1^2-18*p00*x0*x1*y0-18*p00*x0*x1*y1+18*p01*x0*x1*y0+18*p01*x0*x1*y1+18*p10*x0*x1*y0+18*p10*x0*x1*y1-18*p11*x0*x1*y0-18*p11*x0*x1*y1)/((y0-y1)*(x0^3*y0^2-2*x0^3*y0*y1+x0^3*y1^2-3*x0^2*x1*y0^2+6*x0^2*x1*y0*y1-3*x0^2*x1*y1^2+3*x0*x1^2*y0^2-6*x0*x1^2*y0*y1+3*x0*x1^2*y1^2-x1^3*y0^2+2*x1^3*y0*y1-x1^3*y1^2))
		a13 = -(-2*dxy00*x0^2*x1*y0+2*dxy00*x0^2*x1*y1+dxy00*x0*x1^2*y0-dxy00*x0*x1^2*y1+dxy00*x1^3*y0-dxy00*x1^3*y1-2*dxy01*x0^2*x1*y0+2*dxy01*x0^2*x1*y1+dxy01*x0*x1^2*y0-dxy01*x0*x1^2*y1+dxy01*x1^3*y0-dxy01*x1^3*y1-dxy10*x0^3*y0+dxy10*x0^3*y1-dxy10*x0^2*x1*y0+dxy10*x0^2*x1*y1+2*dxy10*x0*x1^2*y0-2*dxy10*x0*x1^2*y1-dxy11*x0^3*y0+dxy11*x0^3*y1-dxy11*x0^2*x1*y0+dxy11*x0^2*x1*y1+2*dxy11*x0*x1^2*y0-2*dxy11*x0*x1^2*y1+4*dx00*x0^2*x1-2*dx00*x0*x1^2-2*dx00*x1^3-4*dx01*x0^2*x1+2*dx01*x0*x1^2+2*dx01*x1^3+2*dx10*x0^3+2*dx10*x0^2*x1-4*dx10*x0*x1^2-2*dx11*x0^3-2*dx11*x0^2*x1+4*dx11*x0*x1^2+6*dy00*x0*x1*y0-6*dy00*x0*x1*y1+6*dy01*x0*x1*y0-6*dy01*x0*x1*y1-6*dy10*x0*x1*y0+6*dy10*x0*x1*y1-6*dy11*x0*x1*y0+6*dy11*x0*x1*y1-12*p00*x0*x1+12*p01*x0*x1+12*p10*x0*x1-12*p11*x0*x1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
		a20 = -(-dxy00*x0^2*y0^2*y1^2+dxy00*x0^2*y0*y1^3-dxy00*x0*x1*y0^2*y1^2+dxy00*x0*x1*y0*y1^3+2*dxy00*x1^2*y0^2*y1^2-2*dxy00*x1^2*y0*y1^3-dxy01*x0^2*y0^3*y1+dxy01*x0^2*y0^2*y1^2-dxy01*x0*x1*y0^3*y1+dxy01*x0*x1*y0^2*y1^2+2*dxy01*x1^2*y0^3*y1-2*dxy01*x1^2*y0^2*y1^2-2*dxy10*x0^2*y0^2*y1^2+2*dxy10*x0^2*y0*y1^3+dxy10*x0*x1*y0^2*y1^2-dxy10*x0*x1*y0*y1^3+dxy10*x1^2*y0^2*y1^2-dxy10*x1^2*y0*y1^3-2*dxy11*x0^2*y0^3*y1+2*dxy11*x0^2*y0^2*y1^2+dxy11*x0*x1*y0^3*y1-dxy11*x0*x1*y0^2*y1^2+dxy11*x1^2*y0^3*y1-dxy11*x1^2*y0^2*y1^2+3*dx00*x0^2*y0*y1^2-dx00*x0^2*y1^3+3*dx00*x0*x1*y0*y1^2-dx00*x0*x1*y1^3-6*dx00*x1^2*y0*y1^2+2*dx00*x1^2*y1^3+dx01*x0^2*y0^3-3*dx01*x0^2*y0^2*y1+dx01*x0*x1*y0^3-3*dx01*x0*x1*y0^2*y1-2*dx01*x1^2*y0^3+6*dx01*x1^2*y0^2*y1+6*dx10*x0^2*y0*y1^2-2*dx10*x0^2*y1^3-3*dx10*x0*x1*y0*y1^2+dx10*x0*x1*y1^3-3*dx10*x1^2*y0*y1^2+dx10*x1^2*y1^3+2*dx11*x0^2*y0^3-6*dx11*x0^2*y0^2*y1-dx11*x0*x1*y0^3+3*dx11*x0*x1*y0^2*y1-dx11*x1^2*y0^3+3*dx11*x1^2*y0^2*y1+3*dy00*x0*y0^2*y1^2-3*dy00*x0*y0*y1^3+3*dy00*x1*y0^2*y1^2-3*dy00*x1*y0*y1^3+3*dy01*x0*y0^3*y1-3*dy01*x0*y0^2*y1^2+3*dy01*x1*y0^3*y1-3*dy01*x1*y0^2*y1^2-3*dy10*x0*y0^2*y1^2+3*dy10*x0*y0*y1^3-3*dy10*x1*y0^2*y1^2+3*dy10*x1*y0*y1^3-3*dy11*x0*y0^3*y1+3*dy11*x0*y0^2*y1^2-3*dy11*x1*y0^3*y1+3*dy11*x1*y0^2*y1^2-9*p00*x0*y0*y1^2+3*p00*x0*y1^3-9*p00*x1*y0*y1^2+3*p00*x1*y1^3-3*p01*x0*y0^3+9*p01*x0*y0^2*y1-3*p01*x1*y0^3+9*p01*x1*y0^2*y1+9*p10*x0*y0*y1^2-3*p10*x0*y1^3+9*p10*x1*y0*y1^2-3*p10*x1*y1^3+3*p11*x0*y0^3-9*p11*x0*y0^2*y1+3*p11*x1*y0^3-9*p11*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
		a21 = (-2*dxy00*x0^2*y0^2*y1+dxy00*x0^2*y0*y1^2+dxy00*x0^2*y1^3-2*dxy00*x0*x1*y0^2*y1+dxy00*x0*x1*y0*y1^2+dxy00*x0*x1*y1^3+4*dxy00*x1^2*y0^2*y1-2*dxy00*x1^2*y0*y1^2-2*dxy00*x1^2*y1^3-dxy01*x0^2*y0^3-dxy01*x0^2*y0^2*y1+2*dxy01*x0^2*y0*y1^2-dxy01*x0*x1*y0^3-dxy01*x0*x1*y0^2*y1+2*dxy01*x0*x1*y0*y1^2+2*dxy01*x1^2*y0^3+2*dxy01*x1^2*y0^2*y1-4*dxy01*x1^2*y0*y1^2-4*dxy10*x0^2*y0^2*y1+2*dxy10*x0^2*y0*y1^2+2*dxy10*x0^2*y1^3+2*dxy10*x0*x1*y0^2*y1-dxy10*x0*x1*y0*y1^2-dxy10*x0*x1*y1^3+2*dxy10*x1^2*y0^2*y1-dxy10*x1^2*y0*y1^2-dxy10*x1^2*y1^3-2*dxy11*x0^2*y0^3-2*dxy11*x0^2*y0^2*y1+4*dxy11*x0^2*y0*y1^2+dxy11*x0*x1*y0^3+dxy11*x0*x1*y0^2*y1-2*dxy11*x0*x1*y0*y1^2+dxy11*x1^2*y0^3+dxy11*x1^2*y0^2*y1-2*dxy11*x1^2*y0*y1^2+6*dx00*x0^2*y0*y1+6*dx00*x0*x1*y0*y1-12*dx00*x1^2*y0*y1-6*dx01*x0^2*y0*y1-6*dx01*x0*x1*y0*y1+12*dx01*x1^2*y0*y1+12*dx10*x0^2*y0*y1-6*dx10*x0*x1*y0*y1-6*dx10*x1^2*y0*y1-12*dx11*x0^2*y0*y1+6*dx11*x0*x1*y0*y1+6*dx11*x1^2*y0*y1+6*dy00*x0*y0^2*y1-3*dy00*x0*y0*y1^2-3*dy00*x0*y1^3+6*dy00*x1*y0^2*y1-3*dy00*x1*y0*y1^2-3*dy00*x1*y1^3+3*dy01*x0*y0^3+3*dy01*x0*y0^2*y1-6*dy01*x0*y0*y1^2+3*dy01*x1*y0^3+3*dy01*x1*y0^2*y1-6*dy01*x1*y0*y1^2-6*dy10*x0*y0^2*y1+3*dy10*x0*y0*y1^2+3*dy10*x0*y1^3-6*dy10*x1*y0^2*y1+3*dy10*x1*y0*y1^2+3*dy10*x1*y1^3-3*dy11*x0*y0^3-3*dy11*x0*y0^2*y1+6*dy11*x0*y0*y1^2-3*dy11*x1*y0^3-3*dy11*x1*y0^2*y1+6*dy11*x1*y0*y1^2-18*p00*x0*y0*y1-18*p00*x1*y0*y1+18*p01*x0*y0*y1+18*p01*x1*y0*y1+18*p10*x0*y0*y1+18*p10*x1*y0*y1-18*p11*x0*y0*y1-18*p11*x1*y0*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
		a22 = -(-dxy00*x0^2*y0^2-dxy00*x0^2*y0*y1+2*dxy00*x0^2*y1^2-dxy00*x0*x1*y0^2-dxy00*x0*x1*y0*y1+2*dxy00*x0*x1*y1^2+2*dxy00*x1^2*y0^2+2*dxy00*x1^2*y0*y1-4*dxy00*x1^2*y1^2-2*dxy01*x0^2*y0^2+dxy01*x0^2*y0*y1+dxy01*x0^2*y1^2-2*dxy01*x0*x1*y0^2+dxy01*x0*x1*y0*y1+dxy01*x0*x1*y1^2+4*dxy01*x1^2*y0^2-2*dxy01*x1^2*y0*y1-2*dxy01*x1^2*y1^2-2*dxy10*x0^2*y0^2-2*dxy10*x0^2*y0*y1+4*dxy10*x0^2*y1^2+dxy10*x0*x1*y0^2+dxy10*x0*x1*y0*y1-2*dxy10*x0*x1*y1^2+dxy10*x1^2*y0^2+dxy10*x1^2*y0*y1-2*dxy10*x1^2*y1^2-4*dxy11*x0^2*y0^2+2*dxy11*x0^2*y0*y1+2*dxy11*x0^2*y1^2+2*dxy11*x0*x1*y0^2-dxy11*x0*x1*y0*y1-dxy11*x0*x1*y1^2+2*dxy11*x1^2*y0^2-dxy11*x1^2*y0*y1-dxy11*x1^2*y1^2+3*dx00*x0^2*y0+3*dx00*x0^2*y1+3*dx00*x0*x1*y0+3*dx00*x0*x1*y1-6*dx00*x1^2*y0-6*dx00*x1^2*y1-3*dx01*x0^2*y0-3*dx01*x0^2*y1-3*dx01*x0*x1*y0-3*dx01*x0*x1*y1+6*dx01*x1^2*y0+6*dx01*x1^2*y1+6*dx10*x0^2*y0+6*dx10*x0^2*y1-3*dx10*x0*x1*y0-3*dx10*x0*x1*y1-3*dx10*x1^2*y0-3*dx10*x1^2*y1-6*dx11*x0^2*y0-6*dx11*x0^2*y1+3*dx11*x0*x1*y0+3*dx11*x0*x1*y1+3*dx11*x1^2*y0+3*dx11*x1^2*y1+3*dy00*x0*y0^2+3*dy00*x0*y0*y1-6*dy00*x0*y1^2+3*dy00*x1*y0^2+3*dy00*x1*y0*y1-6*dy00*x1*y1^2+6*dy01*x0*y0^2-3*dy01*x0*y0*y1-3*dy01*x0*y1^2+6*dy01*x1*y0^2-3*dy01*x1*y0*y1-3*dy01*x1*y1^2-3*dy10*x0*y0^2-3*dy10*x0*y0*y1+6*dy10*x0*y1^2-3*dy10*x1*y0^2-3*dy10*x1*y0*y1+6*dy10*x1*y1^2-6*dy11*x0*y0^2+3*dy11*x0*y0*y1+3*dy11*x0*y1^2-6*dy11*x1*y0^2+3*dy11*x1*y0*y1+3*dy11*x1*y1^2-9*p00*x0*y0-9*p00*x0*y1-9*p00*x1*y0-9*p00*x1*y1+9*p01*x0*y0+9*p01*x0*y1+9*p01*x1*y0+9*p01*x1*y1+9*p10*x0*y0+9*p10*x0*y1+9*p10*x1*y0+9*p10*x1*y1-9*p11*x0*y0-9*p11*x0*y1-9*p11*x1*y0-9*p11*x1*y1)/((x0*y0-x0*y1-x1*y0+x1*y1)*(x0^2*y0^2-2*x0^2*y0*y1+x0^2*y1^2-2*x0*x1*y0^2+4*x0*x1*y0*y1-2*x0*x1*y1^2+x1^2*y0^2-2*x1^2*y0*y1+x1^2*y1^2))
		a23 = (-dxy00*x0^2*y0+dxy00*x0^2*y1-dxy00*x0*x1*y0+dxy00*x0*x1*y1+2*dxy00*x1^2*y0-2*dxy00*x1^2*y1-dxy01*x0^2*y0+dxy01*x0^2*y1-dxy01*x0*x1*y0+dxy01*x0*x1*y1+2*dxy01*x1^2*y0-2*dxy01*x1^2*y1-2*dxy10*x0^2*y0+2*dxy10*x0^2*y1+dxy10*x0*x1*y0-dxy10*x0*x1*y1+dxy10*x1^2*y0-dxy10*x1^2*y1-2*dxy11*x0^2*y0+2*dxy11*x0^2*y1+dxy11*x0*x1*y0-dxy11*x0*x1*y1+dxy11*x1^2*y0-dxy11*x1^2*y1+2*dx00*x0^2+2*dx00*x0*x1-4*dx00*x1^2-2*dx01*x0^2-2*dx01*x0*x1+4*dx01*x1^2+4*dx10*x0^2-2*dx10*x0*x1-2*dx10*x1^2-4*dx11*x0^2+2*dx11*x0*x1+2*dx11*x1^2+3*dy00*x0*y0-3*dy00*x0*y1+3*dy00*x1*y0-3*dy00*x1*y1+3*dy01*x0*y0-3*dy01*x0*y1+3*dy01*x1*y0-3*dy01*x1*y1-3*dy10*x0*y0+3*dy10*x0*y1-3*dy10*x1*y0+3*dy10*x1*y1-3*dy11*x0*y0+3*dy11*x0*y1-3*dy11*x1*y0+3*dy11*x1*y1-6*p00*x0-6*p00*x1+6*p01*x0+6*p01*x1+6*p10*x0+6*p10*x1-6*p11*x0-6*p11*x1)/((x0*y0^3-3*x0*y0^2*y1+3*x0*y0*y1^2-x0*y1^3-x1*y0^3+3*x1*y0^2*y1-3*x1*y0*y1^2+x1*y1^3)*(x0^2-2*x0*x1+x1^2))
		a30 = (-dxy00*x0*y0^2*y1^2+dxy00*x0*y0*y1^3+dxy00*x1*y0^2*y1^2-dxy00*x1*y0*y1^3-dxy01*x0*y0^3*y1+dxy01*x0*y0^2*y1^2+dxy01*x1*y0^3*y1-dxy01*x1*y0^2*y1^2-dxy10*x0*y0^2*y1^2+dxy10*x0*y0*y1^3+dxy10*x1*y0^2*y1^2-dxy10*x1*y0*y1^3-dxy11*x0*y0^3*y1+dxy11*x0*y0^2*y1^2+dxy11*x1*y0^3*y1-dxy11*x1*y0^2*y1^2+3*dx00*x0*y0*y1^2-dx00*x0*y1^3-3*dx00*x1*y0*y1^2+dx00*x1*y1^3+dx01*x0*y0^3-3*dx01*x0*y0^2*y1-dx01*x1*y0^3+3*dx01*x1*y0^2*y1+3*dx10*x0*y0*y1^2-dx10*x0*y1^3-3*dx10*x1*y0*y1^2+dx10*x1*y1^3+dx11*x0*y0^3-3*dx11*x0*y0^2*y1-dx11*x1*y0^3+3*dx11*x1*y0^2*y1+2*dy00*y0^2*y1^2-2*dy00*y0*y1^3+2*dy01*y0^3*y1-2*dy01*y0^2*y1^2-2*dy10*y0^2*y1^2+2*dy10*y0*y1^3-2*dy11*y0^3*y1+2*dy11*y0^2*y1^2-6*p00*y0*y1^2+2*p00*y1^3-2*p01*y0^3+6*p01*y0^2*y1+6*p10*y0*y1^2-2*p10*y1^3+2*p11*y0^3-6*p11*y0^2*y1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
		a31 = -(-2*dxy00*x0*y0^2*y1+dxy00*x0*y0*y1^2+dxy00*x0*y1^3+2*dxy00*x1*y0^2*y1-dxy00*x1*y0*y1^2-dxy00*x1*y1^3-dxy01*x0*y0^3-dxy01*x0*y0^2*y1+2*dxy01*x0*y0*y1^2+dxy01*x1*y0^3+dxy01*x1*y0^2*y1-2*dxy01*x1*y0*y1^2-2*dxy10*x0*y0^2*y1+dxy10*x0*y0*y1^2+dxy10*x0*y1^3+2*dxy10*x1*y0^2*y1-dxy10*x1*y0*y1^2-dxy10*x1*y1^3-dxy11*x0*y0^3-dxy11*x0*y0^2*y1+2*dxy11*x0*y0*y1^2+dxy11*x1*y0^3+dxy11*x1*y0^2*y1-2*dxy11*x1*y0*y1^2+6*dx00*x0*y0*y1-6*dx00*x1*y0*y1-6*dx01*x0*y0*y1+6*dx01*x1*y0*y1+6*dx10*x0*y0*y1-6*dx10*x1*y0*y1-6*dx11*x0*y0*y1+6*dx11*x1*y0*y1+4*dy00*y0^2*y1-2*dy00*y0*y1^2-2*dy00*y1^3+2*dy01*y0^3+2*dy01*y0^2*y1-4*dy01*y0*y1^2-4*dy10*y0^2*y1+2*dy10*y0*y1^2+2*dy10*y1^3-2*dy11*y0^3-2*dy11*y0^2*y1+4*dy11*y0*y1^2-12*p00*y0*y1+12*p01*y0*y1+12*p10*y0*y1-12*p11*y0*y1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
		a32 = (-dxy00*x0*y0^2-dxy00*x0*y0*y1+2*dxy00*x0*y1^2+dxy00*x1*y0^2+dxy00*x1*y0*y1-2*dxy00*x1*y1^2-2*dxy01*x0*y0^2+dxy01*x0*y0*y1+dxy01*x0*y1^2+2*dxy01*x1*y0^2-dxy01*x1*y0*y1-dxy01*x1*y1^2-dxy10*x0*y0^2-dxy10*x0*y0*y1+2*dxy10*x0*y1^2+dxy10*x1*y0^2+dxy10*x1*y0*y1-2*dxy10*x1*y1^2-2*dxy11*x0*y0^2+dxy11*x0*y0*y1+dxy11*x0*y1^2+2*dxy11*x1*y0^2-dxy11*x1*y0*y1-dxy11*x1*y1^2+3*dx00*x0*y0+3*dx00*x0*y1-3*dx00*x1*y0-3*dx00*x1*y1-3*dx01*x0*y0-3*dx01*x0*y1+3*dx01*x1*y0+3*dx01*x1*y1+3*dx10*x0*y0+3*dx10*x0*y1-3*dx10*x1*y0-3*dx10*x1*y1-3*dx11*x0*y0-3*dx11*x0*y1+3*dx11*x1*y0+3*dx11*x1*y1+2*dy00*y0^2+2*dy00*y0*y1-4*dy00*y1^2+4*dy01*y0^2-2*dy01*y0*y1-2*dy01*y1^2-2*dy10*y0^2-2*dy10*y0*y1+4*dy10*y1^2-4*dy11*y0^2+2*dy11*y0*y1+2*dy11*y1^2-6*p00*y0-6*p00*y1+6*p01*y0+6*p01*y1+6*p10*y0+6*p10*y1-6*p11*y0-6*p11*y1)/((x0^3*y0-x0^3*y1-3*x0^2*x1*y0+3*x0^2*x1*y1+3*x0*x1^2*y0-3*x0*x1^2*y1-x1^3*y0+x1^3*y1)*(y0^2-2*y0*y1+y1^2))
		a33 = -(-dxy00*x0*y0+dxy00*x0*y1+dxy00*x1*y0-dxy00*x1*y1-dxy01*x0*y0+dxy01*x0*y1+dxy01*x1*y0-dxy01*x1*y1-dxy10*x0*y0+dxy10*x0*y1+dxy10*x1*y0-dxy10*x1*y1-dxy11*x0*y0+dxy11*x0*y1+dxy11*x1*y0-dxy11*x1*y1+2*dx00*x0-2*dx00*x1-2*dx01*x0+2*dx01*x1+2*dx10*x0-2*dx10*x1-2*dx11*x0+2*dx11*x1+2*dy00*y0-2*dy00*y1+2*dy01*y0-2*dy01*y1-2*dy10*y0+2*dy10*y1-2*dy11*y0+2*dy11*y1-4*p00+4*p01+4*p10-4*p11)/(x0^3*y0^3-3*x0^3*y0^2*y1+3*x0^3*y0*y1^2-x0^3*y1^3-3*x0^2*x1*y0^3+9*x0^2*x1*y0^2*y1-9*x0^2*x1*y0*y1^2+3*x0^2*x1*y1^3+3*x0*x1^2*y0^3-9*x0*x1^2*y0^2*y1+9*x0*x1^2*y0*y1^2-3*x0*x1^2*y1^3-x1^3*y0^3+3*x1^3*y0^2*y1-3*x1^3*y0*y1^2+x1^3*y1^3)

--
\end{code}

{-
A bicubic equation will have the general form 
f(x,y) = a00 + a01*y + a02*y^2 + a03*y^3 + a10*x + a11*x*y + a12*x*y^2 + a13*x*y^3 + a20*(x^2) + a21*(x^2)*y + a22*(x^2)*y^2 + a23*(x^2)*y^3 + a30*(x^3) + a31*(x^3)*y + a32*(x^3)*y^2 + a33*(x^3)*y^3

Given four points with values 
f (x0,y0) = p00 
f (x1,y0) = p10
f (x0,y1) = p01
f (x1,y1) = p11
=> 
f((x0),(y0)) = a00 + a01*(y0) + a02*(y0)^2 + a03*(y0)^3 + a10*(x0) + a11*(x0)*(y0) + a12*(x0)*(y0)^2 + a13*(x0)*(y0)^3 + a20*((x0)^2) + a21*((x0)^2)*(y0) + a22*((x0)^2)*(y0)^2 + a23*((x0)^2)*(y0)^3 + a30*((x0)^3) + a31*((x0)^3)*(y0) + a32*((x0)^3)*(y0)^2 + a33*((x0)^3)*(y0)^3 = p00 
f((x1),(y0)) = a00 + a01*(y0) + a02*(y0)^2 + a03*(y0)^3 + a10*(x1) + a11*(x1)*(y0) + a12*(x1)*(y0)^2 + a13*(x1)*(y0)^3 + a20*((x1)^2) + a21*((x1)^2)*(y0) + a22*((x1)^2)*(y0)^2 + a23*((x1)^2)*(y0)^3 + a30*((x1)^3) + a31*((x1)^3)*(y0) + a32*((x1)^3)*(y0)^2 + a33*((x1)^3)*(y0)^3 = p10 
f((x0),(y1)) = a00 + a01*(y1) + a02*(y1)^2 + a03*(y1)^3 + a10*(x0) + a11*(x0)*(y1) + a12*(x0)*(y1)^2 + a13*(x0)*(y1)^3 + a20*((x0)^2) + a21*((x0)^2)*(y1) + a22*((x0)^2)*(y1)^2 + a23*((x0)^2)*(y1)^3 + a30*((x0)^3) + a31*((x0)^3)*(y1) + a32*((x0)^3)*(y1)^2 + a33*((x0)^3)*(y1)^3 = p01 
f((x1),(y1)) = a00 + a01*(y1) + a02*(y1)^2 + a03*(y1)^3 + a10*(x1) + a11*(x1)*(y1) + a12*(x1)*(y1)^2 + a13*(x1)*(y1)^3 + a20*((x1)^2) + a21*((x1)^2)*(y1) + a22*((x1)^2)*(y1)^2 + a23*((x1)^2)*(y1)^3 + a30*((x1)^3) + a31*((x1)^3)*(y1) + a32*((x1)^3)*(y1)^2 + a33*((x1)^3)*(y1)^3 = p11 


derivatives in x
f'x (x0,y0) = dx00
f'x (x1,y0) = dx10
f'x (x0,y1) = dx01
f'x (x1,y1) = dx11
&
f'x(x,y) = a10 + a11*y + a12*y^2 + a13*y^3 + a20*(2*x) + a21*(2*x)*y + a22*(2*x)*y^2 + a23*(2*x)*y^3 + a30*(3*x^2) + a31*(3*x^2)*y + a32*(3*x^2)*y^2 + a33*(3*x^2)*y^3
=>
f'x((x0),(y0)) = a10 + a11*(y0) + a12*(y0)^2 + a13*(y0)^3 + a20*(2*(x0)) + a21*(2*(x0))*(y0) + a22*(2*(x0))*(y0)^2 + a23*(2*(x0))*(y0)^3 + a30*(3*(x0)^2) + a31*(3*(x0)^2)*(y0) + a32*(3*(x0)^2)*(y0)^2 + a33*(3*(x0)^2)*(y0)^3 = dx00
f'x((x1),(y0)) = a10 + a11*(y0) + a12*(y0)^2 + a13*(y0)^3 + a20*(2*(x1)) + a21*(2*(x1))*(y0) + a22*(2*(x1))*(y0)^2 + a23*(2*(x1))*(y0)^3 + a30*(3*(x1)^2) + a31*(3*(x1)^2)*(y0) + a32*(3*(x1)^2)*(y0)^2 + a33*(3*(x1)^2)*(y0)^3 = dx10
f'x((x0),(y1)) = a10 + a11*(y1) + a12*(y1)^2 + a13*(y1)^3 + a20*(2*(x0)) + a21*(2*(x0))*(y1) + a22*(2*(x0))*(y1)^2 + a23*(2*(x0))*(y1)^3 + a30*(3*(x0)^2) + a31*(3*(x0)^2)*(y1) + a32*(3*(x0)^2)*(y1)^2 + a33*(3*(x0)^2)*(y1)^3 = dx01
f'x((x1),(y1)) = a10 + a11*(y1) + a12*(y1)^2 + a13*(y1)^3 + a20*(2*(x1)) + a21*(2*(x1))*(y1) + a22*(2*(x1))*(y1)^2 + a23*(2*(x1))*(y1)^3 + a30*(3*(x1)^2) + a31*(3*(x1)^2)*(y1) + a32*(3*(x1)^2)*(y1)^2 + a33*(3*(x1)^2)*(y1)^3 = dx11


derivatives in y
f'y (x0,y0) = dy00
f'y (x1,y0) = dy10
f'y (x0,y1) = dy01
f'y (x1,y1) = dy11
&
f'y (x,y) = a01 + a02*2*y + a03*3*y^2 + a11*x + a12*x*2*y + a13*x*3*y^2 + a21*(x^2) + a22*(x^2)*2*y + a23*(x^2)*3*y^2 + a31*(x^3) + a32*(x^3)*2*y + a33*(x^3)*3*y^2
=>
f'y ((x0),(y0)) = a01 + a02*2*(y0) + a03*3*(y0)^2 + a11*(x0) + a12*(x0)*2*(y0) + a13*(x0)*3*(y0)^2 + a21*((x0)^2) + a22*((x0)^2)*2*(y0) + a23*((x0)^2)*3*(y0)^2 + a31*((x0)^3) + a32*((x0)^3)*2*(y0) + a33*((x0)^3)*3*(y0)^2 = dy00
f'y ((x1),(y0)) = a01 + a02*2*(y0) + a03*3*(y0)^2 + a11*(x1) + a12*(x1)*2*(y0) + a13*(x1)*3*(y0)^2 + a21*((x1)^2) + a22*((x1)^2)*2*(y0) + a23*((x1)^2)*3*(y0)^2 + a31*((x1)^3) + a32*((x1)^3)*2*(y0) + a33*((x1)^3)*3*(y0)^2 = dy10
f'y ((x0),(y1)) = a01 + a02*2*(y1) + a03*3*(y1)^2 + a11*(x0) + a12*(x0)*2*(y1) + a13*(x0)*3*(y1)^2 + a21*((x0)^2) + a22*((x0)^2)*2*(y1) + a23*((x0)^2)*3*(y1)^2 + a31*((x0)^3) + a32*((x0)^3)*2*(y1) + a33*((x0)^3)*3*(y1)^2 = dy01
f'y ((x1),(y1)) = a01 + a02*2*(y1) + a03*3*(y1)^2 + a11*(x1) + a12*(x1)*2*(y1) + a13*(x1)*3*(y1)^2 + a21*((x1)^2) + a22*((x1)^2)*2*(y1) + a23*((x1)^2)*3*(y1)^2 + a31*((x1)^3) + a32*((x1)^3)*2*(y1) + a33*((x1)^3)*3*(y1)^2 = dy11


cross derivatives
f''xy (x0,y0) = dxy00
f''xy (x1,y0) = dxy10
f''xy (x0,y1) = dxy01
f''xy (x1,y1) = dxy11
&
f'xy (x,y) = a11 + a12*2*y + a13*3*y^2 + a21*(2*x) + a22*(2*x)*2*y + a23*(2*x)*3*y^2  + a31*(3*x^2) + a32*(3*x^2)*2*y + a33*(3*x^2)*3*y^2
=>
f'xy ((x0),(y0)) = a11 + a12*2*(y0) + a13*3*(y0)^2 + a21*(2*(x0)) + a22*(2*(x0))*2*(y0) + a23*(2*(x0))*3*(y0)^2  + a31*(3*(x0)^2) + a32*(3*(x0)^2)*2*(y0) + a33*(3*(x0)^2)*3*(y0)^2 = dxy00
f'xy ((x1),(y0)) = a11 + a12*2*(y0) + a13*3*(y0)^2 + a21*(2*(x1)) + a22*(2*(x1))*2*(y0) + a23*(2*(x1))*3*(y0)^2  + a31*(3*(x1)^2) + a32*(3*(x1)^2)*2*(y0) + a33*(3*(x1)^2)*3*(y0)^2 = dxy10
f'xy ((x0),(y1)) = a11 + a12*2*(y1) + a13*3*(y1)^2 + a21*(2*(x0)) + a22*(2*(x0))*2*(y1) + a23*(2*(x0))*3*(y1)^2  + a31*(3*(x0)^2) + a32*(3*(x0)^2)*2*(y1) + a33*(3*(x0)^2)*3*(y1)^2 = dxy01
f'xy ((x1),(y1)) = a11 + a12*2*(y1) + a13*3*(y1)^2 + a21*(2*(x1)) + a22*(2*(x1))*2*(y1) + a23*(2*(x1))*3*(y1)^2  + a31*(3*(x1)^2) + a32*(3*(x1)^2)*2*(y1) + a33*(3*(x1)^2)*3*(y1)^2 = dxy11


{a00 + a01*(y0) + a02*(y0)^2 + a03*(y0)^3 + a10*(x0) + a11*(x0)*(y0) + a12*(x0)*(y0)^2 + a13*(x0)*(y0)^3 + a20*((x0)^2) + a21*((x0)^2)*(y0) + a22*((x0)^2)*(y0)^2 + a23*((x0)^2)*(y0)^3 + a30*((x0)^3) + a31*((x0)^3)*(y0) + a32*((x0)^3)*(y0)^2 + a33*((x0)^3)*(y0)^3 = p00,
a00 + a01*(y0) + a02*(y0)^2 + a03*(y0)^3 + a10*(x1) + a11*(x1)*(y0) + a12*(x1)*(y0)^2 + a13*(x1)*(y0)^3 + a20*((x1)^2) + a21*((x1)^2)*(y0) + a22*((x1)^2)*(y0)^2 + a23*((x1)^2)*(y0)^3 + a30*((x1)^3) + a31*((x1)^3)*(y0) + a32*((x1)^3)*(y0)^2 + a33*((x1)^3)*(y0)^3 = p10,
a00 + a01*(y1) + a02*(y1)^2 + a03*(y1)^3 + a10*(x0) + a11*(x0)*(y1) + a12*(x0)*(y1)^2 + a13*(x0)*(y1)^3 + a20*((x0)^2) + a21*((x0)^2)*(y1) + a22*((x0)^2)*(y1)^2 + a23*((x0)^2)*(y1)^3 + a30*((x0)^3) + a31*((x0)^3)*(y1) + a32*((x0)^3)*(y1)^2 + a33*((x0)^3)*(y1)^3 = p01,
a00 + a01*(y1) + a02*(y1)^2 + a03*(y1)^3 + a10*(x1) + a11*(x1)*(y1) + a12*(x1)*(y1)^2 + a13*(x1)*(y1)^3 + a20*((x1)^2) + a21*((x1)^2)*(y1) + a22*((x1)^2)*(y1)^2 + a23*((x1)^2)*(y1)^3 + a30*((x1)^3) + a31*((x1)^3)*(y1) + a32*((x1)^3)*(y1)^2 + a33*((x1)^3)*(y1)^3 = p11,
a10 + a11*(y0) + a12*(y0)^2 + a13*(y0)^3 + a20*(2*(x0)) + a21*(2*(x0))*(y0) + a22*(2*(x0))*(y0)^2 + a23*(2*(x0))*(y0)^3 + a30*(3*(x0)^2) + a31*(3*(x0)^2)*(y0) + a32*(3*(x0)^2)*(y0)^2 + a33*(3*(x0)^2)*(y0)^3 = dx00,
a10 + a11*(y0) + a12*(y0)^2 + a13*(y0)^3 + a20*(2*(x1)) + a21*(2*(x1))*(y0) + a22*(2*(x1))*(y0)^2 + a23*(2*(x1))*(y0)^3 + a30*(3*(x1)^2) + a31*(3*(x1)^2)*(y0) + a32*(3*(x1)^2)*(y0)^2 + a33*(3*(x1)^2)*(y0)^3 = dx10,
a10 + a11*(y1) + a12*(y1)^2 + a13*(y1)^3 + a20*(2*(x0)) + a21*(2*(x0))*(y1) + a22*(2*(x0))*(y1)^2 + a23*(2*(x0))*(y1)^3 + a30*(3*(x0)^2) + a31*(3*(x0)^2)*(y1) + a32*(3*(x0)^2)*(y1)^2 + a33*(3*(x0)^2)*(y1)^3 = dx01,
a10 + a11*(y1) + a12*(y1)^2 + a13*(y1)^3 + a20*(2*(x1)) + a21*(2*(x1))*(y1) + a22*(2*(x1))*(y1)^2 + a23*(2*(x1))*(y1)^3 + a30*(3*(x1)^2) + a31*(3*(x1)^2)*(y1) + a32*(3*(x1)^2)*(y1)^2 + a33*(3*(x1)^2)*(y1)^3 = dx11,
a01 + a02*2*(y0) + a03*3*(y0)^2 + a11*(x0) + a12*(x0)*2*(y0) + a13*(x0)*3*(y0)^2 + a21*((x0)^2) + a22*((x0)^2)*2*(y0) + a23*((x0)^2)*3*(y0)^2 + a31*((x0)^3) + a32*((x0)^3)*2*(y0) + a33*((x0)^3)*3*(y0)^2 = dy00,
a01 + a02*2*(y0) + a03*3*(y0)^2 + a11*(x1) + a12*(x1)*2*(y0) + a13*(x1)*3*(y0)^2 + a21*((x1)^2) + a22*((x1)^2)*2*(y0) + a23*((x1)^2)*3*(y0)^2 + a31*((x1)^3) + a32*((x1)^3)*2*(y0) + a33*((x1)^3)*3*(y0)^2 = dy10,
a01 + a02*2*(y1) + a03*3*(y1)^2 + a11*(x0) + a12*(x0)*2*(y1) + a13*(x0)*3*(y1)^2 + a21*((x0)^2) + a22*((x0)^2)*2*(y1) + a23*((x0)^2)*3*(y1)^2 + a31*((x0)^3) + a32*((x0)^3)*2*(y1) + a33*((x0)^3)*3*(y1)^2 = dy01,
a01 + a02*2*(y1) + a03*3*(y1)^2 + a11*(x1) + a12*(x1)*2*(y1) + a13*(x1)*3*(y1)^2 + a21*((x1)^2) + a22*((x1)^2)*2*(y1) + a23*((x1)^2)*3*(y1)^2 + a31*((x1)^3) + a32*((x1)^3)*2*(y1) + a33*((x1)^3)*3*(y1)^2 = dy11,
a11 + a12*2*(y0) + a13*3*(y0)^2 + a21*(2*(x0)) + a22*(2*(x0))*2*(y0) + a23*(2*(x0))*3*(y0)^2  + a31*(3*(x0)^2) + a32*(3*(x0)^2)*2*(y0) + a33*(3*(x0)^2)*3*(y0)^2 = dxy00,
a11 + a12*2*(y0) + a13*3*(y0)^2 + a21*(2*(x1)) + a22*(2*(x1))*2*(y0) + a23*(2*(x1))*3*(y0)^2  + a31*(3*(x1)^2) + a32*(3*(x1)^2)*2*(y0) + a33*(3*(x1)^2)*3*(y0)^2 = dxy10,
a11 + a12*2*(y1) + a13*3*(y1)^2 + a21*(2*(x0)) + a22*(2*(x0))*2*(y1) + a23*(2*(x0))*3*(y1)^2  + a31*(3*(x0)^2) + a32*(3*(x0)^2)*2*(y1) + a33*(3*(x0)^2)*3*(y1)^2 = dxy01,
a11 + a12*2*(y1) + a13*3*(y1)^2 + a21*(2*(x1)) + a22*(2*(x1))*2*(y1) + a23*(2*(x1))*3*(y1)^2  + a31*(3*(x1)^2) + a32*(3*(x1)^2)*2*(y1) + a33*(3*(x1)^2)*3*(y1)^2 = dxy11}

with params x0,x1,y0,y1,p00,p10,p01,p11,dx00,dx10,dx01,dx11,dy00,dy10,dy01,dy11,dxy00,dxy10,dxy01,dxy11
To be Solved for
{a00,a10,a20,a30 ,a01 ,a11,a21,a31 ,a02 ,a12 ,a22 ,a32 ,a03 ,a13 ,a23 ,a33 }

Giving

a00 = -(-dxy00*x0^2*x1^2*y0^2*y1^2+dxy00*x0^2*x1^2*y0*y1^3+dxy00*x0*x1^3*y0^2*y1^2-dxy00*x0*x1^3*y0*y1^3-dxy01*x0^2*x1^2*y0^3*y1+dxy01*x0^2*x1^2*y0^2*y1^2+dxy01*x0*x1^3*y0^3*y1-dxy01*x0*x1^3*y0^2*y1^2-dxy10*x0^3*x1*y0^2*y1^2+dxy10*x0^3*x1*y0*y1^3+dxy10*x0^2*x1^2*y0^2*y1^2-dxy10*x0^2*x1^2*y0*y1^3-dxy11*x0^3*x1*y0^3*y1+dxy11*x0^3*x1*y0^2*y1^2+dxy11*x0^2*x1^2*y0^3*y1-dxy11*x0^2*x1^2*y0^2*y1^2+3*dx00*x0^2*x1^2*y0*y1^2-dx00*x0^2*x1^2*y1^3-3*dx00*x0*x1^3*y0*y1^2+dx00*x0*x1^3*y1^3+dx01*x0^2*x1^2*y0^3-3*dx01*x0^2*x1^2*y0^2*y1-dx01*x0*x1^3*y0^3+3*dx01*x0*x1^3*y0^2*y1+3*dx10*x0^3*x1*y0*y1^2-dx10*x0^3*x1*y1^3-3*dx10*x0^2*x1^2*y0*y1^2+dx10*x0^2*x1^2*y1^3+dx11*x0^3*x1*y0^3-3*dx11*x0^3*x1*y0^2*y1-dx11*x0^2*x1^2*y0^3+3*dx11*x0^2*x1^2*y0^2*y1+3*dy00*x0*x1^2*y0^2*y1^2-3*dy00*x0*x1^2*y0*y1^3-dy00*x1^3*y0^2*y1^2+dy00*x1^3*y0*y1^3+3*dy01*x0*x1^2*y0^3*y1-3*dy01*x0*x1^2*y0^2*y1^2-dy01*x1^3*y0^3*y1+dy01*x1^3*y0^2*y1^2+dy10*x0^3*y0^2*y1^2-dy10*x0^3*y0*y1^3-3*dy10*x0^2*x1*y0^2*y1^2+3*dy10*x0^2*x1*y0*y1^3+dy11*x0^3*y0^3*y1-dy11*x0^3*y0^2*y1^2-3*dy11*x0^2*x1*y0^3*y1+3*dy11*x0^2*x1*y0^2*y1^2-9*p00*x0*x1^2*y0*y1^2+3*p00*x0*x1^2*y1^3+3*p00*x1^3*y0*y1^2-p00*x1^3*y1^3-3*p01*x0*x1^2*y0^3+9*p01*x0*x1^2*y0^2*y1+p01*x1^3*y0^3-3*p01*x1^3*y0^2*y1-3*p10*x0^3*y0*y1^2+p10*x0^3*y1^3+9*p10*x0^2*x1*y0*y1^2-3*p10*x0^2*x1*y1^3-p11*x0^3*y0^3+3*p11*x0^3*y0^2*y1+3*p11*x0^2*x1*y0^3-9*p11*x0^2*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
a01 = (-2*dxy00*x0^2*x1^2*y0^2*y1+dxy00*x0^2*x1^2*y0*y1^2+dxy00*x0^2*x1^2*y1^3+2*dxy00*x0*x1^3*y0^2*y1-dxy00*x0*x1^3*y0*y1^2-dxy00*x0*x1^3*y1^3-dxy01*x0^2*x1^2*y0^3-dxy01*x0^2*x1^2*y0^2*y1+2*dxy01*x0^2*x1^2*y0*y1^2+dxy01*x0*x1^3*y0^3+dxy01*x0*x1^3*y0^2*y1-2*dxy01*x0*x1^3*y0*y1^2-2*dxy10*x0^3*x1*y0^2*y1+dxy10*x0^3*x1*y0*y1^2+dxy10*x0^3*x1*y1^3+2*dxy10*x0^2*x1^2*y0^2*y1-dxy10*x0^2*x1^2*y0*y1^2-dxy10*x0^2*x1^2*y1^3-dxy11*x0^3*x1*y0^3-dxy11*x0^3*x1*y0^2*y1+2*dxy11*x0^3*x1*y0*y1^2+dxy11*x0^2*x1^2*y0^3+dxy11*x0^2*x1^2*y0^2*y1-2*dxy11*x0^2*x1^2*y0*y1^2+6*dx00*x0^2*x1^2*y0*y1-6*dx00*x0*x1^3*y0*y1-6*dx01*x0^2*x1^2*y0*y1+6*dx01*x0*x1^3*y0*y1+6*dx10*x0^3*x1*y0*y1-6*dx10*x0^2*x1^2*y0*y1-6*dx11*x0^3*x1*y0*y1+6*dx11*x0^2*x1^2*y0*y1+6*dy00*x0*x1^2*y0^2*y1-3*dy00*x0*x1^2*y0*y1^2-3*dy00*x0*x1^2*y1^3-2*dy00*x1^3*y0^2*y1+dy00*x1^3*y0*y1^2+dy00*x1^3*y1^3+3*dy01*x0*x1^2*y0^3+3*dy01*x0*x1^2*y0^2*y1-6*dy01*x0*x1^2*y0*y1^2-dy01*x1^3*y0^3-dy01*x1^3*y0^2*y1+2*dy01*x1^3*y0*y1^2+2*dy10*x0^3*y0^2*y1-dy10*x0^3*y0*y1^2-dy10*x0^3*y1^3-6*dy10*x0^2*x1*y0^2*y1+3*dy10*x0^2*x1*y0*y1^2+3*dy10*x0^2*x1*y1^3+dy11*x0^3*y0^3+dy11*x0^3*y0^2*y1-2*dy11*x0^3*y0*y1^2-3*dy11*x0^2*x1*y0^3-3*dy11*x0^2*x1*y0^2*y1+6*dy11*x0^2*x1*y0*y1^2-18*p00*x0*x1^2*y0*y1+6*p00*x1^3*y0*y1+18*p01*x0*x1^2*y0*y1-6*p01*x1^3*y0*y1-6*p10*x0^3*y0*y1+18*p10*x0^2*x1*y0*y1+6*p11*x0^3*y0*y1-18*p11*x0^2*x1*y0*y1)/((x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3)*(x0-x1))
a02 = -(-dxy00*x0^2*x1^2*y0^2-dxy00*x0^2*x1^2*y0*y1+2*dxy00*x0^2*x1^2*y1^2+dxy00*x0*x1^3*y0^2+dxy00*x0*x1^3*y0*y1-2*dxy00*x0*x1^3*y1^2-2*dxy01*x0^2*x1^2*y0^2+dxy01*x0^2*x1^2*y0*y1+dxy01*x0^2*x1^2*y1^2+2*dxy01*x0*x1^3*y0^2-dxy01*x0*x1^3*y0*y1-dxy01*x0*x1^3*y1^2-dxy10*x0^3*x1*y0^2-dxy10*x0^3*x1*y0*y1+2*dxy10*x0^3*x1*y1^2+dxy10*x0^2*x1^2*y0^2+dxy10*x0^2*x1^2*y0*y1-2*dxy10*x0^2*x1^2*y1^2-2*dxy11*x0^3*x1*y0^2+dxy11*x0^3*x1*y0*y1+dxy11*x0^3*x1*y1^2+2*dxy11*x0^2*x1^2*y0^2-dxy11*x0^2*x1^2*y0*y1-dxy11*x0^2*x1^2*y1^2+3*dx00*x0^2*x1^2*y0+3*dx00*x0^2*x1^2*y1-3*dx00*x0*x1^3*y0-3*dx00*x0*x1^3*y1-3*dx01*x0^2*x1^2*y0-3*dx01*x0^2*x1^2*y1+3*dx01*x0*x1^3*y0+3*dx01*x0*x1^3*y1+3*dx10*x0^3*x1*y0+3*dx10*x0^3*x1*y1-3*dx10*x0^2*x1^2*y0-3*dx10*x0^2*x1^2*y1-3*dx11*x0^3*x1*y0-3*dx11*x0^3*x1*y1+3*dx11*x0^2*x1^2*y0+3*dx11*x0^2*x1^2*y1+3*dy00*x0*x1^2*y0^2+3*dy00*x0*x1^2*y0*y1-6*dy00*x0*x1^2*y1^2-dy00*x1^3*y0^2-dy00*x1^3*y0*y1+2*dy00*x1^3*y1^2+6*dy01*x0*x1^2*y0^2-3*dy01*x0*x1^2*y0*y1-3*dy01*x0*x1^2*y1^2-2*dy01*x1^3*y0^2+dy01*x1^3*y0*y1+dy01*x1^3*y1^2+dy10*x0^3*y0^2+dy10*x0^3*y0*y1-2*dy10*x0^3*y1^2-3*dy10*x0^2*x1*y0^2-3*dy10*x0^2*x1*y0*y1+6*dy10*x0^2*x1*y1^2+2*dy11*x0^3*y0^2-dy11*x0^3*y0*y1-dy11*x0^3*y1^2-6*dy11*x0^2*x1*y0^2+3*dy11*x0^2*x1*y0*y1+3*dy11*x0^2*x1*y1^2-9*p00*x0*x1^2*y0-9*p00*x0*x1^2*y1+3*p00*x1^3*y0+3*p00*x1^3*y1+9*p01*x0*x1^2*y0+9*p01*x0*x1^2*y1-3*p01*x1^3*y0-3*p01*x1^3*y1-3*p10*x0^3*y0-3*p10*x0^3*y1+9*p10*x0^2*x1*y0+9*p10*x0^2*x1*y1+3*p11*x0^3*y0+3*p11*x0^3*y1-9*p11*x0^2*x1*y0-9*p11*x0^2*x1*y1)/((y0-y1)*(x0^3*y0^2-2*x0^3*y0*y1+x0^3*y1^2-3*x0^2*x1*y0^2+6*x0^2*x1*y0*y1-3*x0^2*x1*y1^2+3*x0*x1^2*y0^2-6*x0*x1^2*y0*y1+3*x0*x1^2*y1^2-x1^3*y0^2+2*x1^3*y0*y1-x1^3*y1^2))
a03 = (-dxy00*x0^2*x1^2*y0+dxy00*x0^2*x1^2*y1+dxy00*x0*x1^3*y0-dxy00*x0*x1^3*y1-dxy01*x0^2*x1^2*y0+dxy01*x0^2*x1^2*y1+dxy01*x0*x1^3*y0-dxy01*x0*x1^3*y1-dxy10*x0^3*x1*y0+dxy10*x0^3*x1*y1+dxy10*x0^2*x1^2*y0-dxy10*x0^2*x1^2*y1-dxy11*x0^3*x1*y0+dxy11*x0^3*x1*y1+dxy11*x0^2*x1^2*y0-dxy11*x0^2*x1^2*y1+2*dx00*x0^2*x1^2-2*dx00*x0*x1^3-2*dx01*x0^2*x1^2+2*dx01*x0*x1^3+2*dx10*x0^3*x1-2*dx10*x0^2*x1^2-2*dx11*x0^3*x1+2*dx11*x0^2*x1^2+3*dy00*x0*x1^2*y0-3*dy00*x0*x1^2*y1-dy00*x1^3*y0+dy00*x1^3*y1+3*dy01*x0*x1^2*y0-3*dy01*x0*x1^2*y1-dy01*x1^3*y0+dy01*x1^3*y1+dy10*x0^3*y0-dy10*x0^3*y1-3*dy10*x0^2*x1*y0+3*dy10*x0^2*x1*y1+dy11*x0^3*y0-dy11*x0^3*y1-3*dy11*x0^2*x1*y0+3*dy11*x0^2*x1*y1-6*p00*x0*x1^2+2*p00*x1^3+6*p01*x0*x1^2-2*p01*x1^3-2*p10*x0^3+6*p10*x0^2*x1+2*p11*x0^3-6*p11*x0^2*x1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
a10 = (-2*dxy00*x0^2*x1*y0^2*y1^2+2*dxy00*x0^2*x1*y0*y1^3+dxy00*x0*x1^2*y0^2*y1^2-dxy00*x0*x1^2*y0*y1^3+dxy00*x1^3*y0^2*y1^2-dxy00*x1^3*y0*y1^3-2*dxy01*x0^2*x1*y0^3*y1+2*dxy01*x0^2*x1*y0^2*y1^2+dxy01*x0*x1^2*y0^3*y1-dxy01*x0*x1^2*y0^2*y1^2+dxy01*x1^3*y0^3*y1-dxy01*x1^3*y0^2*y1^2-dxy10*x0^3*y0^2*y1^2+dxy10*x0^3*y0*y1^3-dxy10*x0^2*x1*y0^2*y1^2+dxy10*x0^2*x1*y0*y1^3+2*dxy10*x0*x1^2*y0^2*y1^2-2*dxy10*x0*x1^2*y0*y1^3-dxy11*x0^3*y0^3*y1+dxy11*x0^3*y0^2*y1^2-dxy11*x0^2*x1*y0^3*y1+dxy11*x0^2*x1*y0^2*y1^2+2*dxy11*x0*x1^2*y0^3*y1-2*dxy11*x0*x1^2*y0^2*y1^2+6*dx00*x0^2*x1*y0*y1^2-2*dx00*x0^2*x1*y1^3-3*dx00*x0*x1^2*y0*y1^2+dx00*x0*x1^2*y1^3-3*dx00*x1^3*y0*y1^2+dx00*x1^3*y1^3+2*dx01*x0^2*x1*y0^3-6*dx01*x0^2*x1*y0^2*y1-dx01*x0*x1^2*y0^3+3*dx01*x0*x1^2*y0^2*y1-dx01*x1^3*y0^3+3*dx01*x1^3*y0^2*y1+3*dx10*x0^3*y0*y1^2-dx10*x0^3*y1^3+3*dx10*x0^2*x1*y0*y1^2-dx10*x0^2*x1*y1^3-6*dx10*x0*x1^2*y0*y1^2+2*dx10*x0*x1^2*y1^3+dx11*x0^3*y0^3-3*dx11*x0^3*y0^2*y1+dx11*x0^2*x1*y0^3-3*dx11*x0^2*x1*y0^2*y1-2*dx11*x0*x1^2*y0^3+6*dx11*x0*x1^2*y0^2*y1+6*dy00*x0*x1*y0^2*y1^2-6*dy00*x0*x1*y0*y1^3+6*dy01*x0*x1*y0^3*y1-6*dy01*x0*x1*y0^2*y1^2-6*dy10*x0*x1*y0^2*y1^2+6*dy10*x0*x1*y0*y1^3-6*dy11*x0*x1*y0^3*y1+6*dy11*x0*x1*y0^2*y1^2-18*p00*x0*x1*y0*y1^2+6*p00*x0*x1*y1^3-6*p01*x0*x1*y0^3+18*p01*x0*x1*y0^2*y1+18*p10*x0*x1*y0*y1^2-6*p10*x0*x1*y1^3+6*p11*x0*x1*y0^3-18*p11*x0*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
a11 = -(-4*dxy00*x0^2*x1*y0^2*y1+2*dxy00*x0^2*x1*y0*y1^2+2*dxy00*x0^2*x1*y1^3+2*dxy00*x0*x1^2*y0^2*y1-dxy00*x0*x1^2*y0*y1^2-dxy00*x0*x1^2*y1^3+2*dxy00*x1^3*y0^2*y1-dxy00*x1^3*y0*y1^2-dxy00*x1^3*y1^3-2*dxy01*x0^2*x1*y0^3-2*dxy01*x0^2*x1*y0^2*y1+4*dxy01*x0^2*x1*y0*y1^2+dxy01*x0*x1^2*y0^3+dxy01*x0*x1^2*y0^2*y1-2*dxy01*x0*x1^2*y0*y1^2+dxy01*x1^3*y0^3+dxy01*x1^3*y0^2*y1-2*dxy01*x1^3*y0*y1^2-2*dxy10*x0^3*y0^2*y1+dxy10*x0^3*y0*y1^2+dxy10*x0^3*y1^3-2*dxy10*x0^2*x1*y0^2*y1+dxy10*x0^2*x1*y0*y1^2+dxy10*x0^2*x1*y1^3+4*dxy10*x0*x1^2*y0^2*y1-2*dxy10*x0*x1^2*y0*y1^2-2*dxy10*x0*x1^2*y1^3-dxy11*x0^3*y0^3-dxy11*x0^3*y0^2*y1+2*dxy11*x0^3*y0*y1^2-dxy11*x0^2*x1*y0^3-dxy11*x0^2*x1*y0^2*y1+2*dxy11*x0^2*x1*y0*y1^2+2*dxy11*x0*x1^2*y0^3+2*dxy11*x0*x1^2*y0^2*y1-4*dxy11*x0*x1^2*y0*y1^2+12*dx00*x0^2*x1*y0*y1-6*dx00*x0*x1^2*y0*y1-6*dx00*x1^3*y0*y1-12*dx01*x0^2*x1*y0*y1+6*dx01*x0*x1^2*y0*y1+6*dx01*x1^3*y0*y1+6*dx10*x0^3*y0*y1+6*dx10*x0^2*x1*y0*y1-12*dx10*x0*x1^2*y0*y1-6*dx11*x0^3*y0*y1-6*dx11*x0^2*x1*y0*y1+12*dx11*x0*x1^2*y0*y1+12*dy00*x0*x1*y0^2*y1-6*dy00*x0*x1*y0*y1^2-6*dy00*x0*x1*y1^3+6*dy01*x0*x1*y0^3+6*dy01*x0*x1*y0^2*y1-12*dy01*x0*x1*y0*y1^2-12*dy10*x0*x1*y0^2*y1+6*dy10*x0*x1*y0*y1^2+6*dy10*x0*x1*y1^3-6*dy11*x0*x1*y0^3-6*dy11*x0*x1*y0^2*y1+12*dy11*x0*x1*y0*y1^2-36*p00*x0*x1*y0*y1+36*p01*x0*x1*y0*y1+36*p10*x0*x1*y0*y1-36*p11*x0*x1*y0*y1)/((x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3)*(x0-x1))
a12 = (-2*dxy00*x0^2*x1*y0^2-2*dxy00*x0^2*x1*y0*y1+4*dxy00*x0^2*x1*y1^2+dxy00*x0*x1^2*y0^2+dxy00*x0*x1^2*y0*y1-2*dxy00*x0*x1^2*y1^2+dxy00*x1^3*y0^2+dxy00*x1^3*y0*y1-2*dxy00*x1^3*y1^2-4*dxy01*x0^2*x1*y0^2+2*dxy01*x0^2*x1*y0*y1+2*dxy01*x0^2*x1*y1^2+2*dxy01*x0*x1^2*y0^2-dxy01*x0*x1^2*y0*y1-dxy01*x0*x1^2*y1^2+2*dxy01*x1^3*y0^2-dxy01*x1^3*y0*y1-dxy01*x1^3*y1^2-dxy10*x0^3*y0^2-dxy10*x0^3*y0*y1+2*dxy10*x0^3*y1^2-dxy10*x0^2*x1*y0^2-dxy10*x0^2*x1*y0*y1+2*dxy10*x0^2*x1*y1^2+2*dxy10*x0*x1^2*y0^2+2*dxy10*x0*x1^2*y0*y1-4*dxy10*x0*x1^2*y1^2-2*dxy11*x0^3*y0^2+dxy11*x0^3*y0*y1+dxy11*x0^3*y1^2-2*dxy11*x0^2*x1*y0^2+dxy11*x0^2*x1*y0*y1+dxy11*x0^2*x1*y1^2+4*dxy11*x0*x1^2*y0^2-2*dxy11*x0*x1^2*y0*y1-2*dxy11*x0*x1^2*y1^2+6*dx00*x0^2*x1*y0+6*dx00*x0^2*x1*y1-3*dx00*x0*x1^2*y0-3*dx00*x0*x1^2*y1-3*dx00*x1^3*y0-3*dx00*x1^3*y1-6*dx01*x0^2*x1*y0-6*dx01*x0^2*x1*y1+3*dx01*x0*x1^2*y0+3*dx01*x0*x1^2*y1+3*dx01*x1^3*y0+3*dx01*x1^3*y1+3*dx10*x0^3*y0+3*dx10*x0^3*y1+3*dx10*x0^2*x1*y0+3*dx10*x0^2*x1*y1-6*dx10*x0*x1^2*y0-6*dx10*x0*x1^2*y1-3*dx11*x0^3*y0-3*dx11*x0^3*y1-3*dx11*x0^2*x1*y0-3*dx11*x0^2*x1*y1+6*dx11*x0*x1^2*y0+6*dx11*x0*x1^2*y1+6*dy00*x0*x1*y0^2+6*dy00*x0*x1*y0*y1-12*dy00*x0*x1*y1^2+12*dy01*x0*x1*y0^2-6*dy01*x0*x1*y0*y1-6*dy01*x0*x1*y1^2-6*dy10*x0*x1*y0^2-6*dy10*x0*x1*y0*y1+12*dy10*x0*x1*y1^2-12*dy11*x0*x1*y0^2+6*dy11*x0*x1*y0*y1+6*dy11*x0*x1*y1^2-18*p00*x0*x1*y0-18*p00*x0*x1*y1+18*p01*x0*x1*y0+18*p01*x0*x1*y1+18*p10*x0*x1*y0+18*p10*x0*x1*y1-18*p11*x0*x1*y0-18*p11*x0*x1*y1)/((y0-y1)*(x0^3*y0^2-2*x0^3*y0*y1+x0^3*y1^2-3*x0^2*x1*y0^2+6*x0^2*x1*y0*y1-3*x0^2*x1*y1^2+3*x0*x1^2*y0^2-6*x0*x1^2*y0*y1+3*x0*x1^2*y1^2-x1^3*y0^2+2*x1^3*y0*y1-x1^3*y1^2))
a13 = -(-2*dxy00*x0^2*x1*y0+2*dxy00*x0^2*x1*y1+dxy00*x0*x1^2*y0-dxy00*x0*x1^2*y1+dxy00*x1^3*y0-dxy00*x1^3*y1-2*dxy01*x0^2*x1*y0+2*dxy01*x0^2*x1*y1+dxy01*x0*x1^2*y0-dxy01*x0*x1^2*y1+dxy01*x1^3*y0-dxy01*x1^3*y1-dxy10*x0^3*y0+dxy10*x0^3*y1-dxy10*x0^2*x1*y0+dxy10*x0^2*x1*y1+2*dxy10*x0*x1^2*y0-2*dxy10*x0*x1^2*y1-dxy11*x0^3*y0+dxy11*x0^3*y1-dxy11*x0^2*x1*y0+dxy11*x0^2*x1*y1+2*dxy11*x0*x1^2*y0-2*dxy11*x0*x1^2*y1+4*dx00*x0^2*x1-2*dx00*x0*x1^2-2*dx00*x1^3-4*dx01*x0^2*x1+2*dx01*x0*x1^2+2*dx01*x1^3+2*dx10*x0^3+2*dx10*x0^2*x1-4*dx10*x0*x1^2-2*dx11*x0^3-2*dx11*x0^2*x1+4*dx11*x0*x1^2+6*dy00*x0*x1*y0-6*dy00*x0*x1*y1+6*dy01*x0*x1*y0-6*dy01*x0*x1*y1-6*dy10*x0*x1*y0+6*dy10*x0*x1*y1-6*dy11*x0*x1*y0+6*dy11*x0*x1*y1-12*p00*x0*x1+12*p01*x0*x1+12*p10*x0*x1-12*p11*x0*x1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
a20 = -(-dxy00*x0^2*y0^2*y1^2+dxy00*x0^2*y0*y1^3-dxy00*x0*x1*y0^2*y1^2+dxy00*x0*x1*y0*y1^3+2*dxy00*x1^2*y0^2*y1^2-2*dxy00*x1^2*y0*y1^3-dxy01*x0^2*y0^3*y1+dxy01*x0^2*y0^2*y1^2-dxy01*x0*x1*y0^3*y1+dxy01*x0*x1*y0^2*y1^2+2*dxy01*x1^2*y0^3*y1-2*dxy01*x1^2*y0^2*y1^2-2*dxy10*x0^2*y0^2*y1^2+2*dxy10*x0^2*y0*y1^3+dxy10*x0*x1*y0^2*y1^2-dxy10*x0*x1*y0*y1^3+dxy10*x1^2*y0^2*y1^2-dxy10*x1^2*y0*y1^3-2*dxy11*x0^2*y0^3*y1+2*dxy11*x0^2*y0^2*y1^2+dxy11*x0*x1*y0^3*y1-dxy11*x0*x1*y0^2*y1^2+dxy11*x1^2*y0^3*y1-dxy11*x1^2*y0^2*y1^2+3*dx00*x0^2*y0*y1^2-dx00*x0^2*y1^3+3*dx00*x0*x1*y0*y1^2-dx00*x0*x1*y1^3-6*dx00*x1^2*y0*y1^2+2*dx00*x1^2*y1^3+dx01*x0^2*y0^3-3*dx01*x0^2*y0^2*y1+dx01*x0*x1*y0^3-3*dx01*x0*x1*y0^2*y1-2*dx01*x1^2*y0^3+6*dx01*x1^2*y0^2*y1+6*dx10*x0^2*y0*y1^2-2*dx10*x0^2*y1^3-3*dx10*x0*x1*y0*y1^2+dx10*x0*x1*y1^3-3*dx10*x1^2*y0*y1^2+dx10*x1^2*y1^3+2*dx11*x0^2*y0^3-6*dx11*x0^2*y0^2*y1-dx11*x0*x1*y0^3+3*dx11*x0*x1*y0^2*y1-dx11*x1^2*y0^3+3*dx11*x1^2*y0^2*y1+3*dy00*x0*y0^2*y1^2-3*dy00*x0*y0*y1^3+3*dy00*x1*y0^2*y1^2-3*dy00*x1*y0*y1^3+3*dy01*x0*y0^3*y1-3*dy01*x0*y0^2*y1^2+3*dy01*x1*y0^3*y1-3*dy01*x1*y0^2*y1^2-3*dy10*x0*y0^2*y1^2+3*dy10*x0*y0*y1^3-3*dy10*x1*y0^2*y1^2+3*dy10*x1*y0*y1^3-3*dy11*x0*y0^3*y1+3*dy11*x0*y0^2*y1^2-3*dy11*x1*y0^3*y1+3*dy11*x1*y0^2*y1^2-9*p00*x0*y0*y1^2+3*p00*x0*y1^3-9*p00*x1*y0*y1^2+3*p00*x1*y1^3-3*p01*x0*y0^3+9*p01*x0*y0^2*y1-3*p01*x1*y0^3+9*p01*x1*y0^2*y1+9*p10*x0*y0*y1^2-3*p10*x0*y1^3+9*p10*x1*y0*y1^2-3*p10*x1*y1^3+3*p11*x0*y0^3-9*p11*x0*y0^2*y1+3*p11*x1*y0^3-9*p11*x1*y0^2*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
a21 = (-2*dxy00*x0^2*y0^2*y1+dxy00*x0^2*y0*y1^2+dxy00*x0^2*y1^3-2*dxy00*x0*x1*y0^2*y1+dxy00*x0*x1*y0*y1^2+dxy00*x0*x1*y1^3+4*dxy00*x1^2*y0^2*y1-2*dxy00*x1^2*y0*y1^2-2*dxy00*x1^2*y1^3-dxy01*x0^2*y0^3-dxy01*x0^2*y0^2*y1+2*dxy01*x0^2*y0*y1^2-dxy01*x0*x1*y0^3-dxy01*x0*x1*y0^2*y1+2*dxy01*x0*x1*y0*y1^2+2*dxy01*x1^2*y0^3+2*dxy01*x1^2*y0^2*y1-4*dxy01*x1^2*y0*y1^2-4*dxy10*x0^2*y0^2*y1+2*dxy10*x0^2*y0*y1^2+2*dxy10*x0^2*y1^3+2*dxy10*x0*x1*y0^2*y1-dxy10*x0*x1*y0*y1^2-dxy10*x0*x1*y1^3+2*dxy10*x1^2*y0^2*y1-dxy10*x1^2*y0*y1^2-dxy10*x1^2*y1^3-2*dxy11*x0^2*y0^3-2*dxy11*x0^2*y0^2*y1+4*dxy11*x0^2*y0*y1^2+dxy11*x0*x1*y0^3+dxy11*x0*x1*y0^2*y1-2*dxy11*x0*x1*y0*y1^2+dxy11*x1^2*y0^3+dxy11*x1^2*y0^2*y1-2*dxy11*x1^2*y0*y1^2+6*dx00*x0^2*y0*y1+6*dx00*x0*x1*y0*y1-12*dx00*x1^2*y0*y1-6*dx01*x0^2*y0*y1-6*dx01*x0*x1*y0*y1+12*dx01*x1^2*y0*y1+12*dx10*x0^2*y0*y1-6*dx10*x0*x1*y0*y1-6*dx10*x1^2*y0*y1-12*dx11*x0^2*y0*y1+6*dx11*x0*x1*y0*y1+6*dx11*x1^2*y0*y1+6*dy00*x0*y0^2*y1-3*dy00*x0*y0*y1^2-3*dy00*x0*y1^3+6*dy00*x1*y0^2*y1-3*dy00*x1*y0*y1^2-3*dy00*x1*y1^3+3*dy01*x0*y0^3+3*dy01*x0*y0^2*y1-6*dy01*x0*y0*y1^2+3*dy01*x1*y0^3+3*dy01*x1*y0^2*y1-6*dy01*x1*y0*y1^2-6*dy10*x0*y0^2*y1+3*dy10*x0*y0*y1^2+3*dy10*x0*y1^3-6*dy10*x1*y0^2*y1+3*dy10*x1*y0*y1^2+3*dy10*x1*y1^3-3*dy11*x0*y0^3-3*dy11*x0*y0^2*y1+6*dy11*x0*y0*y1^2-3*dy11*x1*y0^3-3*dy11*x1*y0^2*y1+6*dy11*x1*y0*y1^2-18*p00*x0*y0*y1-18*p00*x1*y0*y1+18*p01*x0*y0*y1+18*p01*x1*y0*y1+18*p10*x0*y0*y1+18*p10*x1*y0*y1-18*p11*x0*y0*y1-18*p11*x1*y0*y1)/((x0-x1)*(x0^2*y0^3-3*x0^2*y0^2*y1+3*x0^2*y0*y1^2-x0^2*y1^3-2*x0*x1*y0^3+6*x0*x1*y0^2*y1-6*x0*x1*y0*y1^2+2*x0*x1*y1^3+x1^2*y0^3-3*x1^2*y0^2*y1+3*x1^2*y0*y1^2-x1^2*y1^3))
a22 = -(-dxy00*x0^2*y0^2-dxy00*x0^2*y0*y1+2*dxy00*x0^2*y1^2-dxy00*x0*x1*y0^2-dxy00*x0*x1*y0*y1+2*dxy00*x0*x1*y1^2+2*dxy00*x1^2*y0^2+2*dxy00*x1^2*y0*y1-4*dxy00*x1^2*y1^2-2*dxy01*x0^2*y0^2+dxy01*x0^2*y0*y1+dxy01*x0^2*y1^2-2*dxy01*x0*x1*y0^2+dxy01*x0*x1*y0*y1+dxy01*x0*x1*y1^2+4*dxy01*x1^2*y0^2-2*dxy01*x1^2*y0*y1-2*dxy01*x1^2*y1^2-2*dxy10*x0^2*y0^2-2*dxy10*x0^2*y0*y1+4*dxy10*x0^2*y1^2+dxy10*x0*x1*y0^2+dxy10*x0*x1*y0*y1-2*dxy10*x0*x1*y1^2+dxy10*x1^2*y0^2+dxy10*x1^2*y0*y1-2*dxy10*x1^2*y1^2-4*dxy11*x0^2*y0^2+2*dxy11*x0^2*y0*y1+2*dxy11*x0^2*y1^2+2*dxy11*x0*x1*y0^2-dxy11*x0*x1*y0*y1-dxy11*x0*x1*y1^2+2*dxy11*x1^2*y0^2-dxy11*x1^2*y0*y1-dxy11*x1^2*y1^2+3*dx00*x0^2*y0+3*dx00*x0^2*y1+3*dx00*x0*x1*y0+3*dx00*x0*x1*y1-6*dx00*x1^2*y0-6*dx00*x1^2*y1-3*dx01*x0^2*y0-3*dx01*x0^2*y1-3*dx01*x0*x1*y0-3*dx01*x0*x1*y1+6*dx01*x1^2*y0+6*dx01*x1^2*y1+6*dx10*x0^2*y0+6*dx10*x0^2*y1-3*dx10*x0*x1*y0-3*dx10*x0*x1*y1-3*dx10*x1^2*y0-3*dx10*x1^2*y1-6*dx11*x0^2*y0-6*dx11*x0^2*y1+3*dx11*x0*x1*y0+3*dx11*x0*x1*y1+3*dx11*x1^2*y0+3*dx11*x1^2*y1+3*dy00*x0*y0^2+3*dy00*x0*y0*y1-6*dy00*x0*y1^2+3*dy00*x1*y0^2+3*dy00*x1*y0*y1-6*dy00*x1*y1^2+6*dy01*x0*y0^2-3*dy01*x0*y0*y1-3*dy01*x0*y1^2+6*dy01*x1*y0^2-3*dy01*x1*y0*y1-3*dy01*x1*y1^2-3*dy10*x0*y0^2-3*dy10*x0*y0*y1+6*dy10*x0*y1^2-3*dy10*x1*y0^2-3*dy10*x1*y0*y1+6*dy10*x1*y1^2-6*dy11*x0*y0^2+3*dy11*x0*y0*y1+3*dy11*x0*y1^2-6*dy11*x1*y0^2+3*dy11*x1*y0*y1+3*dy11*x1*y1^2-9*p00*x0*y0-9*p00*x0*y1-9*p00*x1*y0-9*p00*x1*y1+9*p01*x0*y0+9*p01*x0*y1+9*p01*x1*y0+9*p01*x1*y1+9*p10*x0*y0+9*p10*x0*y1+9*p10*x1*y0+9*p10*x1*y1-9*p11*x0*y0-9*p11*x0*y1-9*p11*x1*y0-9*p11*x1*y1)/((x0*y0-x0*y1-x1*y0+x1*y1)*(x0^2*y0^2-2*x0^2*y0*y1+x0^2*y1^2-2*x0*x1*y0^2+4*x0*x1*y0*y1-2*x0*x1*y1^2+x1^2*y0^2-2*x1^2*y0*y1+x1^2*y1^2))
a23 = (-dxy00*x0^2*y0+dxy00*x0^2*y1-dxy00*x0*x1*y0+dxy00*x0*x1*y1+2*dxy00*x1^2*y0-2*dxy00*x1^2*y1-dxy01*x0^2*y0+dxy01*x0^2*y1-dxy01*x0*x1*y0+dxy01*x0*x1*y1+2*dxy01*x1^2*y0-2*dxy01*x1^2*y1-2*dxy10*x0^2*y0+2*dxy10*x0^2*y1+dxy10*x0*x1*y0-dxy10*x0*x1*y1+dxy10*x1^2*y0-dxy10*x1^2*y1-2*dxy11*x0^2*y0+2*dxy11*x0^2*y1+dxy11*x0*x1*y0-dxy11*x0*x1*y1+dxy11*x1^2*y0-dxy11*x1^2*y1+2*dx00*x0^2+2*dx00*x0*x1-4*dx00*x1^2-2*dx01*x0^2-2*dx01*x0*x1+4*dx01*x1^2+4*dx10*x0^2-2*dx10*x0*x1-2*dx10*x1^2-4*dx11*x0^2+2*dx11*x0*x1+2*dx11*x1^2+3*dy00*x0*y0-3*dy00*x0*y1+3*dy00*x1*y0-3*dy00*x1*y1+3*dy01*x0*y0-3*dy01*x0*y1+3*dy01*x1*y0-3*dy01*x1*y1-3*dy10*x0*y0+3*dy10*x0*y1-3*dy10*x1*y0+3*dy10*x1*y1-3*dy11*x0*y0+3*dy11*x0*y1-3*dy11*x1*y0+3*dy11*x1*y1-6*p00*x0-6*p00*x1+6*p01*x0+6*p01*x1+6*p10*x0+6*p10*x1-6*p11*x0-6*p11*x1)/((x0*y0^3-3*x0*y0^2*y1+3*x0*y0*y1^2-x0*y1^3-x1*y0^3+3*x1*y0^2*y1-3*x1*y0*y1^2+x1*y1^3)*(x0^2-2*x0*x1+x1^2))
a30 = (-dxy00*x0*y0^2*y1^2+dxy00*x0*y0*y1^3+dxy00*x1*y0^2*y1^2-dxy00*x1*y0*y1^3-dxy01*x0*y0^3*y1+dxy01*x0*y0^2*y1^2+dxy01*x1*y0^3*y1-dxy01*x1*y0^2*y1^2-dxy10*x0*y0^2*y1^2+dxy10*x0*y0*y1^3+dxy10*x1*y0^2*y1^2-dxy10*x1*y0*y1^3-dxy11*x0*y0^3*y1+dxy11*x0*y0^2*y1^2+dxy11*x1*y0^3*y1-dxy11*x1*y0^2*y1^2+3*dx00*x0*y0*y1^2-dx00*x0*y1^3-3*dx00*x1*y0*y1^2+dx00*x1*y1^3+dx01*x0*y0^3-3*dx01*x0*y0^2*y1-dx01*x1*y0^3+3*dx01*x1*y0^2*y1+3*dx10*x0*y0*y1^2-dx10*x0*y1^3-3*dx10*x1*y0*y1^2+dx10*x1*y1^3+dx11*x0*y0^3-3*dx11*x0*y0^2*y1-dx11*x1*y0^3+3*dx11*x1*y0^2*y1+2*dy00*y0^2*y1^2-2*dy00*y0*y1^3+2*dy01*y0^3*y1-2*dy01*y0^2*y1^2-2*dy10*y0^2*y1^2+2*dy10*y0*y1^3-2*dy11*y0^3*y1+2*dy11*y0^2*y1^2-6*p00*y0*y1^2+2*p00*y1^3-2*p01*y0^3+6*p01*y0^2*y1+6*p10*y0*y1^2-2*p10*y1^3+2*p11*y0^3-6*p11*y0^2*y1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
a31 = -(-2*dxy00*x0*y0^2*y1+dxy00*x0*y0*y1^2+dxy00*x0*y1^3+2*dxy00*x1*y0^2*y1-dxy00*x1*y0*y1^2-dxy00*x1*y1^3-dxy01*x0*y0^3-dxy01*x0*y0^2*y1+2*dxy01*x0*y0*y1^2+dxy01*x1*y0^3+dxy01*x1*y0^2*y1-2*dxy01*x1*y0*y1^2-2*dxy10*x0*y0^2*y1+dxy10*x0*y0*y1^2+dxy10*x0*y1^3+2*dxy10*x1*y0^2*y1-dxy10*x1*y0*y1^2-dxy10*x1*y1^3-dxy11*x0*y0^3-dxy11*x0*y0^2*y1+2*dxy11*x0*y0*y1^2+dxy11*x1*y0^3+dxy11*x1*y0^2*y1-2*dxy11*x1*y0*y1^2+6*dx00*x0*y0*y1-6*dx00*x1*y0*y1-6*dx01*x0*y0*y1+6*dx01*x1*y0*y1+6*dx10*x0*y0*y1-6*dx10*x1*y0*y1-6*dx11*x0*y0*y1+6*dx11*x1*y0*y1+4*dy00*y0^2*y1-2*dy00*y0*y1^2-2*dy00*y1^3+2*dy01*y0^3+2*dy01*y0^2*y1-4*dy01*y0*y1^2-4*dy10*y0^2*y1+2*dy10*y0*y1^2+2*dy10*y1^3-2*dy11*y0^3-2*dy11*y0^2*y1+4*dy11*y0*y1^2-12*p00*y0*y1+12*p01*y0*y1+12*p10*y0*y1-12*p11*y0*y1)/((x0^3-3*x0^2*x1+3*x0*x1^2-x1^3)*(y0^3-3*y0^2*y1+3*y0*y1^2-y1^3))
a32 = (-dxy00*x0*y0^2-dxy00*x0*y0*y1+2*dxy00*x0*y1^2+dxy00*x1*y0^2+dxy00*x1*y0*y1-2*dxy00*x1*y1^2-2*dxy01*x0*y0^2+dxy01*x0*y0*y1+dxy01*x0*y1^2+2*dxy01*x1*y0^2-dxy01*x1*y0*y1-dxy01*x1*y1^2-dxy10*x0*y0^2-dxy10*x0*y0*y1+2*dxy10*x0*y1^2+dxy10*x1*y0^2+dxy10*x1*y0*y1-2*dxy10*x1*y1^2-2*dxy11*x0*y0^2+dxy11*x0*y0*y1+dxy11*x0*y1^2+2*dxy11*x1*y0^2-dxy11*x1*y0*y1-dxy11*x1*y1^2+3*dx00*x0*y0+3*dx00*x0*y1-3*dx00*x1*y0-3*dx00*x1*y1-3*dx01*x0*y0-3*dx01*x0*y1+3*dx01*x1*y0+3*dx01*x1*y1+3*dx10*x0*y0+3*dx10*x0*y1-3*dx10*x1*y0-3*dx10*x1*y1-3*dx11*x0*y0-3*dx11*x0*y1+3*dx11*x1*y0+3*dx11*x1*y1+2*dy00*y0^2+2*dy00*y0*y1-4*dy00*y1^2+4*dy01*y0^2-2*dy01*y0*y1-2*dy01*y1^2-2*dy10*y0^2-2*dy10*y0*y1+4*dy10*y1^2-4*dy11*y0^2+2*dy11*y0*y1+2*dy11*y1^2-6*p00*y0-6*p00*y1+6*p01*y0+6*p01*y1+6*p10*y0+6*p10*y1-6*p11*y0-6*p11*y1)/((x0^3*y0-x0^3*y1-3*x0^2*x1*y0+3*x0^2*x1*y1+3*x0*x1^2*y0-3*x0*x1^2*y1-x1^3*y0+x1^3*y1)*(y0^2-2*y0*y1+y1^2))
a33 = -(-dxy00*x0*y0+dxy00*x0*y1+dxy00*x1*y0-dxy00*x1*y1-dxy01*x0*y0+dxy01*x0*y1+dxy01*x1*y0-dxy01*x1*y1-dxy10*x0*y0+dxy10*x0*y1+dxy10*x1*y0-dxy10*x1*y1-dxy11*x0*y0+dxy11*x0*y1+dxy11*x1*y0-dxy11*x1*y1+2*dx00*x0-2*dx00*x1-2*dx01*x0+2*dx01*x1+2*dx10*x0-2*dx10*x1-2*dx11*x0+2*dx11*x1+2*dy00*y0-2*dy00*y1+2*dy01*y0-2*dy01*y1-2*dy10*y0+2*dy10*y1-2*dy11*y0+2*dy11*y1-4*p00+4*p01+4*p10-4*p11)/(x0^3*y0^3-3*x0^3*y0^2*y1+3*x0^3*y0*y1^2-x0^3*y1^3-3*x0^2*x1*y0^3+9*x0^2*x1*y0^2*y1-9*x0^2*x1*y0*y1^2+3*x0^2*x1*y1^3+3*x0*x1^2*y0^3-9*x0*x1^2*y0^2*y1+9*x0*x1^2*y0*y1^2-3*x0*x1^2*y1^3-x1^3*y0^3+3*x1^3*y0^2*y1-3*x1^3*y0*y1^2+x1^3*y1^3)




-}




2 Dimensional Splines, done analagously to 1d. 

\begin{code}



\end{code}



























super hacky array constructer, for when we just want to plug in values
--topleft (y,x) indexing at 1
\begin{code}

listToArray21 :: [[a]] -> (Int,Int) -> Array (Int,Int) a 
listToArray21 list (n,m) = array ((1,1),(n,m)) [ ((i,j), ( list !! (i - 1) ) !! (j - 1)  ) | i <- [1..n], j <- [1..m]  ]

\end{code}

--bottomleft (x,y) indexing at 0
\begin{code}
listToArray20 :: [[a]] -> (Int,Int) -> Array (Int,Int) a 
listToArray20 list (n,m) = array ((0,0),(n,m)) [ ((i,j), ( (reverse list) !! j ) !! i  ) | i <- [0..n], j <- [0..m]  ]


test1 = listToArray20 [[0,0,0],[0,1.0,0]] (2,1)
test2 = listToArray21 [[0,0,0],[0,1.0,0]] (2,3)

test3 = listToArray20
	[[0,0,0,0,0]
	,[0,0,1.0,0,0]
	,[0,2,3,2,0]
	,[0,-1,1,3,0]
	,[0,0,0,0,0]] (4,4)
	
	
test4 = listToArray20
	[[0,0,0,0,0]
	,[0,0,1.0,0,0]
	,[0,0,0,0,0]
	,[0,0,0,0,0]
	,[0,0,0,0,0]] (4,4)

testfn u = evalPoly2 ( biLinear test3 u ) u

testfn2 u = evalPoly2 ( biCubic test3 u ) u


\end{code}







