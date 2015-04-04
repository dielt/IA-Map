\begin{code}

module Map where


import Util.Base

import Data.Vector as V hiding (empty,elem,map,foldl',foldr)
import Data.Array.Repa hiding (map,foldl',foldr)
import Data.Array.Repa.IO.BMP
import Data.HashMap.Strict hiding (map,foldl',null,foldr)
import Data.List (nub,foldl')


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
	NorthWest -> ["nw","Nw","nW","Northwest","NorthWest","northWest","northwest"]
	NorthEast -> []
	SouthWest -> []
	SouthEast -> []
	Up -> []
	Down -> []
	In -> []
	Out -> []

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

flatArrayToMap :: V.Vector ( V.Vector (Room,[Direction])) -> Map
flatArrayToMap arr = if V.null arr then empty else
	foldl' f empty [0..n]
		where
			getRoom (u,v) = ( arr V.! (fromIntegral u) ) V.! (fromIntegral v)
			n = fromIntegral $ (V.length arr) - 1
			f hmp' x' = foldl' (g m x') hmp' [0..m]
				where m = fromIntegral $ (V.length (arr V.! (fromIntegral x'))) - 1
			g m x hmp y = insert (name r) (addSynonyms (map (h m x y) dirs) r) hmp
				where (r,dirs) = getRoom (x,y)
			h m x y dir = Token ( name . fst . getRoom $ addTuple (x,y) (cardinal2Vec dir) , directionStrings dir)

\end{code}