\begin{code}


module Util.Interpolation where

import Util.Base

import qualified Math.Polynomial as P
import Control.Monad
import Data.List


\end{code}


1 dimensional interpolation, done with tuples

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




Cubic-hermite, finite differance tangents.
\begin{code}

appendPoints :: Fractional a => [(a,a)] -> [(a,a)]
appendPoints xs@(x:y:ys) = 
	let 
		u = last xs
		v = last $ init xs
	in [ ( (fst x) - 1 , linearIntrp x y ((fst x) - 1) ) ] 
		++ xs 
		++ [( (fst u) + 1 , linearIntrp u v ((fst u) + 1)) ] 
appendPoints _ = [] 

--this is not working as intended
hermite :: (Ord a , Fractional a) => [(a,a)] -> Spline a
hermite points = foldr
		(\n s  -> ( \t -> (x n <= t) && (t <= x (n+1)) , P.composePoly ( P.scalePoly (1 / (x (n+1) - x n)) (P.addPoly P.x $ P.constPoly (x n))) ( P.sumPolys 
			[ P.scalePoly (p n) h00 
			, P.scalePoly ((x (n+1) - x n )*(m n)) h10
			, P.scalePoly (p (n+1)) h01
			, P.scalePoly ((x (n+1) - x n )*(m (n+1))) h11
			] )   ) : s 
			) [] [1..((length points) - 1)]
		where
			points' = appendPoints points
			x k = fst $ points' !! k
			p k = snd $ points' !! k
			m k = ( p (k+1) - p k ) / (2*(x (k+1) - x k )) + ( p k - p (k-1) ) / (2*(x k - x (k -1) )) 


--the hermite basis polynomials.
h00 :: (Num a,Eq a) => P.Poly a
h00 = P.sumPolys [ P.scalePoly 2 $ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly (-3) $ P.multPoly P.x P.x , P.one ]
h10 :: (Num a,Eq a) => P.Poly a
h10 = P.sumPolys [ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly (-2) $ P.multPoly P.x P.x , P.x ]
h01 :: (Num a,Eq a) => P.Poly a
h01 = P.sumPolys [ P.scalePoly (-2) $ P.multPoly P.x $ P.multPoly P.x P.x , P.scalePoly (3) $ P.multPoly P.x P.x  ]
h11 :: (Num a,Eq a) => P.Poly a
h11 = P.sumPolys [ P.multPoly P.x $ P.multPoly P.x P.x , P.negatePoly $ P.multPoly P.x P.x  ]

\end{code}
