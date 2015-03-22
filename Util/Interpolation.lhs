\begin{code}


module Util.Interpolation where

import Util.Base

import qualified Math.Polynomial as P
import Control.Monad


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

--
--linear fit by least square
--linearFit :: (Eq a , Fractional a) => [(a,a)] -> P.Poly a
linearFit :: (Eq a , Fractional a) => [(a,a)] -> P.Poly a
linearFit points = 
	let 
		(x',y') = centerFit points
		b = (sum . map (\(a,b) -> (a - x')*(b - y') ) $ points) / (sum . map (\a -> (a - x')^2 ) . map snd $ points )
		a = y' - b*x'
	in P.addPoly (P.scalePoly b P.x) (P.constPoly a)
	

centerFit :: Fractional a => [(a,a)] -> (a,a)
centerFit points = (  (sum (map fst points)) / ( fromIntegral $ length points) , (sum (map snd points)) / (fromIntegral $ length points) )
--
\end{code}



n-th degree polynomials, lagrange method
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










