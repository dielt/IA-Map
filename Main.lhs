
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
main = do
	putStr ">"
	hFlush stdout
	x <- getLine
	let y = read x :: Integer
	let z = testArray3 y
	writeFile "output.pgm" $ arrayToPBM z
	
	
	
\end{code}



