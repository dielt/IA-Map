\begin{code}

module Map where


import Util.Base
import Util.Interpolation

import qualified Data.Array.Unboxed as V 
import qualified Data.Array.Repa as R 
import Data.Array.Repa.IO.BMP
import Data.HashMap.Strict hiding (map,foldl',null,foldr)
import Data.List (nub,foldl')
import qualified Data.Vector.Unboxed as U 
import Data.Maybe


\end{code}

\begin{code}
newtype Token a b = Token (a,[b]) deriving (Eq,Show)
\end{code}

This whole definition is very tentitive, 
The reason why we have a list of functions and not just of strings is that, for example
we will want "north" to translate into differant things in each case.
\begin{code}

data Room = Room 
	{ name :: String
	, adjacent :: [String]
	, synomyms :: [(String -> Maybe String)]
	}

emptyRoom :: String -> Room
emptyRoom s = Room
	{ name = s
	, adjacent = []
	, synomyms = []
	}

addAdjacent :: String -> Room -> Room
addAdjacent str r = r{adjacent = str : (adjacent r)}

addSynonyms' :: Token String String -> Room -> Room
addSynonyms' (Token (adj,syn)) r = 
	r{synomyms = 
		(\x -> if x `elem` syn then Just adj else Nothing
		) : (synomyms r) }

nubAdjacent r  = r{adjacent = nub (adjacent r)}

addSynonym :: Token String String -> Room -> Room
addSynonym t@(Token (adj,_)) = nubAdjacent . addAdjacent adj . addSynonyms' t

addSynonyms :: [Token String String] -> Room -> Room
addSynonyms l z = foldr addSynonym z l



\end{code}


Again tentative
\begin{code}

type Map = HashMap String Room

addRoom :: Map -> Room -> Map
addRoom m r = insert (name r) r m

\end{code}


some basic direction types
\begin{code}

data Direction = North 
	| South 
	| East 
	| West
	| NorthWest
	| NorthEast
	| SouthWest
	| SouthEast
	| Up
	| Down
	| In
	| Out deriving (Eq,Show)
--

directionStrings :: Direction -> [String]
directionStrings d = case d of
	North -> ["n","N","north","North"]
	South -> ["s","S","south","South"]
	East  -> ["e","E","east","East"]
	West -> ["w","W","west","West"]
	NorthWest -> ["nw","Nw","nW","NW","Northwest","NorthWest","northWest","northwest"]
	NorthEast -> ["ne","Ne","nE","NE","Northeast","NorthEast","northEast","northest"]
	SouthWest -> ["sw","Sw","sW","SW","Sorthwest","SorthWest","sorthWest","sorthwest"]
	SouthEast -> ["se","Se","sE","SE","Sortheast","SorthEast","sorthEast","sorthest"]
	Up -> ["u","U","up","Up"]
	Down -> ["d","D","down","Down"]
	In -> ["in","In"]
	Out -> ["out","Out"]

--

cardinal2Vec :: Num a => Direction -> (a,a)
cardinal2Vec d = case d of
	North -> (0,1)
	South -> (0,-1)
	East  -> (1,0)
	West -> (-1,0)
	NorthWest -> (-1,1)
	NorthEast -> (1,1)
	SouthWest -> (-1,-1)
	SouthEast -> (1,-1)
	otherwise -> undefined



\end{code}


\begin{code}
--I'm not sure why this needs to be Integer instead of Int
flatArrayToMap :: V.Array (Integer,Integer) (Room,[Direction]) -> Map
flatArrayToMap arr = foldl' f empty [n_..n'] 
		where
			getRoom (u,v) = arr V.! ((fromIntegral u),(fromIntegral v))
			((n_,n'),(m_,m')) = V.bounds arr
			f hmp' x' = foldl' (g x') hmp' [m_..m']
			g x hmp y = insert (name r) (addSynonyms (map (h x y) dirs) r) hmp
				where (r,dirs) = getRoom (x,y)
			h x y dir = Token ( name . fst . getRoom $ addTuple (x,y) (cardinal2Vec dir) , directionStrings dir)
--


\end{code}


\begin{code}

arrayMap :: (V.Ix a,Enum a, Num a) =>  (b -> c) -> V.Array (a,a) b -> V.Array (a,a) c
arrayMap f arr = V.array ((n_,m_),(n',m')) [((x,y),f (arr V.! (x,y))) | x <- [n_..(n')], y <- [m_..(m')] ]
	where	
		((n_,m_),(n',m')) = V.bounds arr

seedToRandArray ::  Integer -> ((Int,Int),(Int,Int)) -> Int -> V.Array (Int,Int) Int
seedToRandArray seed ((x_,y_),(x',y')) h' = V.array ((x_,y_),(x',y')) [ ((x,y),simpleHashList [x,y,fromIntegral seed] h' ) | x <-[x_..(x')], y <-[y_..(y')] ]

fromIntegralArray :: (V.Ix c,Num c,Enum c,Num b,Integral a) => V.Array (c,c) a -> V.Array (c,c) b
fromIntegralArray = arrayMap fromIntegral


--note we assume n_ and m_ are both 0
embiggenArray :: V.Array (Int,Int) Int -> Int ->  V.Array (Int,Int) Int
embiggenArray arr scale =  V.array ((n_,m_),(n',m')) [ ((x,y) , f (x,y)  ) |  x <- [n_..n'], y <- [m_..m'] ]
	where	
		((n_,m_),(n'',m'')) = V.bounds arr
		n' = (n'')*scale + scale
		m' = (m'')*scale + scale
		splne = bigBicubicSpline $ fromIntegralArray arr
		f (u,v) = round . fromJust $ evalSpline2 splne ((fromIntegral u) / (fromIntegral scale)  , ( fromIntegral v) / (fromIntegral scale)  )
--
testArray :: Integer -> V.Array (Int,Int) Double
testArray seed = fromIntegralArray $ seedToRandArray seed ((0,0),(10,10)) 100 

--testArray2 :: Integer -> V.Array (Int,Int) Int
testArray2 seed = embiggenArray ( seedToRandArray seed ((0,0),(10,10)) 100 ) 10

arrayToRepa :: (U.Unbox a) => V.Array (Int,Int) a -> R.Array R.U R.DIM2 a
arrayToRepa arr = R.fromUnboxed (R.Z R.:. (n' -  1 ) R.:.(m' -  1)) (U.fromList list)
	where	
		((n_,m_),(n',m')) = V.bounds arr
		list = foldr (\a b -> foldr (\x y -> (arr V.! (a,x)) : y ) b [0..m'] ) [] [0..n']




\end{code}


