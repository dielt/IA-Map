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

Quadratic
Continuous first derivatives
\begin{code}

quadraticSpline :: (Ord a , Fractional a) => [(a,a)] -> Spline a
quadraticSpline (x:y:ys) = snd $ foldl' f (y,[(containmentTest x y , linear x y)]) ys
	where
		--f::((a,a),Spline a) -> (a,a) -> ((a,a),Spline a)
		f (u,s) v = ( v ,(containmentTest u v , p ) : s)
			where --see below -- this could be cleaned up
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
				
				
quadraticSpline _ = []



\end{code}





