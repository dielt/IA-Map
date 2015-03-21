\begin{code}


module Util.Interpolation where


import qualified Math.Polynomial as P


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

--linear fit by least square
linearFit :: (Eq a , Fractional a) => [(a,a)] -> P.Poly a
linearFit points = 
	let 
		(x',y') = centerFit points
		b = (sum . map (\(a,b) -> (a - x')*(b - y') ) $ points) / (sum . map (\a -> (a - x')^2 ) . map snd $ points )
		a = y' - b*x'
	in P.addPoly (P.scalePoly b P.x) (P.constPoly a)
	

centerFit :: Fractional a => [(a,a)] -> (a,a)
centerFit points = (  (sum (map fst points)) / ( fromIntegral $ length points) , (sum (map snd points)) / (fromIntegral $ length points) )

\end{code}



n-th degree polynomials, lagrange on n points
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

\end{code}

--spline interpolation, quadratic 
\begin{code}


\end{code}



