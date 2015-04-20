
\begin{code}
module Main where


--import System.Console.ANSI --this would only be used for text bolding


\end{code}


\begin{code}


import Map

import System.Environment
import System.IO

import Data.Array.Repa ()
import Data.Array.Repa.IO.BMP
import Data.Array.Repa.IO.Matrix

\end{code}


\begin{code}
main :: IO()
main = do {-
	putStr ">"
	hFlush stdout
	x <- getLine
	let y = read x :: Integer
	let z1 = test1 y
	let z2 = test2 y
	writeFile "output1.pgm" $ arrayToPBM z1
	writeFile "output2.pgm" $ arrayToPBM z2 --}
	writeFile "output0.pgm" . arrayToPBM  $ test8 9
	writeFile "output1.pgm" . arrayToPBM  $ test9 $ test8 9
	writeFile "output2.pgm" . arrayToPBM  $ test10 $ test8 9
	writeFile "output3.ppm" . arrayToPPM $ heightToColour $ test8 10
	writeFile "output4.ppm" . arrayToPPM $ heightToColour $ test8 11
	writeFile "output5.ppm" . arrayToPPM $ heightToColour $ test8 12
	writeFile "output6.ppm" . arrayToPPM $ heightToColour $ test8 13
	writeFile "output7.ppm" . arrayToPPM $ heightToColour $ test8 14
	
	
	
	
\end{code}



