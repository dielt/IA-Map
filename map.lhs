\begin{code}

module Map where


import Util.Base
import Util.Interpolation

import qualified Data.Array.Unboxed as V 
import qualified Data.Array.Repa as R 
import Data.Array.Repa.IO.BMP
import Data.HashMap.Strict hiding (map,foldl',null,foldr,filter)
import Data.List (nub,foldl')
import qualified Data.Vector.Unboxed as U 
import Data.Maybe
import Data.Monoid


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

arrayMap :: (V.Ix a,Enum a) =>  (b -> c) -> V.Array (a,a) b -> V.Array (a,a) c
arrayMap f arr = V.array ((n_,m_),(n',m')) [((x,y),f (arr V.! (x,y))) | x <- [n_..(n')], y <- [m_..(m')] ]
	where	
		((n_,m_),(n',m')) = V.bounds arr

seedToRandArray :: (V.Ix a,Enum a,Integral a,Integral b) => b -> ((a,a),(a,a)) -> b -> V.Array (a,a) b
seedToRandArray seed ((x_,y_),(x',y')) h' = V.array ((x_,y_),(x',y')) [ ((x,y), simpleHashList [seed,seed,fromIntegral x,fromIntegral y,seed] h' ) | x <-[x_..x'], y <-[y_..y'] ]

fromIntegralArray :: (V.Ix c,Num c,Enum c,Num b,Integral a) => V.Array (c,c) a -> V.Array (c,c) b
fromIntegralArray = arrayMap fromIntegral


clampArray (a,b) = arrayMap (\x -> if x < a then a else if x > b then b else x)

--note we assume n_ and m_ are both 0
embiggenArray :: (V.Ix a,Enum a,Integral a,Integral b) =>  b -> V.Array (a,a) b ->  V.Array (a,a) b
embiggenArray scale arr =  V.array ((n_,m_),(n',m')) [ ((x,y) , f (x,y)  ) |  x <- [n_..n'], y <- [m_..m'] ]
	where	
		((n_,m_),(n'',m'')) = V.bounds arr
		n' = (n'')*(fromIntegral scale)
		m' = (m'')*(fromIntegral scale)
		splne = bigBicubicSpline $ fromIntegralArray arr
		f (u,v) = round . fromJust $ evalSpline2 splne ((fromIntegral u) / (fromIntegral scale)  , ( fromIntegral v) / (fromIntegral scale)  )
--

arrayToRepa :: (V.Ix b,Enum b,Integral b,U.Unbox a) => V.Array (b,b) a -> R.Array R.U R.DIM2 a
arrayToRepa arr = R.fromUnboxed (R.Z R.:. (fromIntegral $ n' -  1 ) R.:.(fromIntegral $ m' -  1)) (U.fromList list)
	where	
		((n_,m_),(n',m')) = V.bounds arr
		list = foldr (\a b -> foldr (\x y -> (arr V.! (a,x)) : y ) b [0..m'] ) [] [0..n']

arrayToList :: (V.Ix b,Enum b) => V.Array (b,b) a -> [[a]]
arrayToList arr = foldr f [] [m_..m']
	where
		((n_,m_),(n',m')) = V.bounds arr
		f m l1 = (foldr (g m) [] [n_..n']) : l1
		g m n l2 = (arr V.! (n,m)) : l2
--

arrayToString :: (V.Ix b,Enum b,Num b,Show a) => V.Array (b,b) a -> String
arrayToString =  unlines . map unwords . arrayToList . arrayMap show

arrayToPBM :: (V.Ix b,Enum b,Num b,Show b,Show a,Num a,Ord a) => V.Array (b,b) a -> String
arrayToPBM arr = "P2\n" ++ dim ++ "\n" ++ height ++ "\n" ++ arrList
	where
		((n_,m_),(n',m')) = V.bounds arr
		dim = show (n' - n_ + 1) ++ " " ++ show (m' - m_ + 1)
		height = show . maximum . concat . arrayToList $ arr
		arrList = arrayToString arr
--

trimArray1 :: (V.Ix b,Enum b,Num a) => V.Array (b,b) a -> V.Array (b,b) a
trimArray1 arr = ((((arr V.// [((x,m_),0)| x <- [n_..n']  ]) V.// [((x,m'),0)| x <- [n_..n']  ]) V.// [((n_,y),0)| y <- [m_..m']  ]) V.// [((n',y),0)| y <- [m_..m']  ])
	where
		((n_,m_),(n',m')) = V.bounds arr

trimArray2 :: (V.Ix b,Enum b,Integral b,Num a) => V.Array (b,b) a -> V.Array (b,b) a
trimArray2 arr = arr V.//  [((x,y),0) | x <- [n_..n'] , y <- [m_..m'] ,  ( ((fromIntegral x) - xc)^2 / a^2 ) + ( ((fromIntegral y) - yc)^2 / b^2  ) > 1  ]
	where
		((n_,m_),(n',m')) = V.bounds arr
		xc = (fromIntegral n_) + a
		yc = (fromIntegral m_) + b
		a = ( fromIntegral (n' - n_) ) / 2
		b = ( fromIntegral (m' - m_) ) / 2
--
addArrays :: (V.Ix a,Enum a,Num b) =>  V.Array (a,a) b -> V.Array (a,a) b  -> V.Array (a,a) b 
addArrays arr1 arr2 = if ((n_,m_),(n',m')) /= V.bounds arr2 then undefined else
	V.array ((n_,m_),(n',m')) [((x,y), (arr1 V.! (x,y)) + (arr2 V.! (x,y)) ) | x <- [n_..n'] , y <- [m_..m'] ]
		where
			((n_,m_),(n',m')) = V.bounds arr1
--

constArray :: (V.Ix a,Enum a) => ((a,a),(a,a)) -> c -> V.Array (a,a) c
constArray ((n_,m_),(n',m')) c = V.array ((n_,m_),(n',m')) [((x,y), c ) | x <- [n_..n'] , y <- [m_..m'] ]

--

sumArrays :: (V.Ix a,Enum a,Num a,Num b) => a -> [V.Array (a,a) b] -> V.Array (a,a) b 
sumArrays dim = foldr addArrays (constArray ((0,0),(dim,dim)) 0 )


test3 = constArray ((0,0),(10,10)) 1

test4' = embiggenArray 2 $ constArray ((0,0),(5,5)) 1

--I really need a better word for this,
fractalInterp :: (V.Ix a,Enum a,Integral a,Integral b) => b ->  [(a,(b,b))] -> a -> [V.Array (a,a) b]
fractalInterp seed ls dim' = 
	map (\(dim_,(h_,h')) -> 
		if dim' `rem` dim_ /= 0 
			then undefined
			else let scale = fromIntegral $ dim' `quot` dim_ in 
				embiggenArray scale . arrayMap ((+) h_) $ seedToRandArray seed ((0,0),(dim_,dim_)) (h' - h_)
	) ls
--
test1 seed = fractalInterp seed [(6,(0,100))] 120

--test2 seed = (clampArray (0,100)) . (arrayMap ((+) (-25))) $  fractalInterp seed [(5,(0,60)),(10,(0,30)),(20,(0,10))] 120

--}
--
normalizeArrayInt :: (V.Ix a,Enum a,Integral b) => b -> V.Array (a,a) b -> V.Array (a,a) b
normalizeArrayInt h arr = arrayMap ( round . (\x -> ( fromIntegral $ x * h)  / m) ) $ arr
	where m = fromIntegral . maximum . concat . arrayToList $ arr :: Double
--}
normalizeArray :: (V.Ix a,Enum a,RealFrac b) => b -> V.Array (a,a) b -> V.Array (a,a) b
normalizeArray h arr = arrayMap ( (\x -> ( x * h)  / m) ) $ arr
	where m = maximum . concat . arrayToList $ arr

--

test4 size = (arrayMap ((+) (1))) . trimArray2 . trimArray1 $ constArray ((0,0),(size,size)) 1

test5 seed = trimArray2 . trimArray1  $ seedToRandArray seed ((0,0),(3,3)) 20

test6 seed = trimArray2 . trimArray1 $ seedToRandArray seed ((0,0),(9,9)) 10

test7 seed = trimArray2 . trimArray1 $ seedToRandArray seed ((0,0),(15,15)) 5


test8 seed = trimArray1 . clampArray (0,15) . arrayMap ((+) (-5)) . normalizeArrayInt 20 $ sumArrays 90 [arr1,arr2,arr3]
	where
		arr1 = embiggenArray 30 $ test5 seed
		arr2 = embiggenArray 10 $ test6 ((seed + 3)^5 + seed^2)
		arr3 = embiggenArray 6 $ test7 (seed^3)
--


\end{code}


\begin{code}

floodFill :: (V.Ix a,Enum a,Integral a) => (b -> Bool) -> (a,a) -> V.Array (a,a) b -> V.Array (a,a) Bool
floodFill test x arr  = floodFillHelper test x arr (( constArray ((n_,m_),(n',m')) False ) V.// [(x,True)] ) []
	where 
		((n_,m_),(n',m')) = V.bounds arr
--

floodFillHelper :: (V.Ix a,Enum a,Ord a,Integral a) => (b -> Bool) -> (a,a) -> V.Array (a,a) b -> V.Array (a,a) Bool -> [(a,a)] -> V.Array (a,a) Bool
floodFillHelper test x arr boolArr stack = 
	let 
		newtest u = and [(test $ arr V.! u) , (not $ boolArr V.! u) ]
		adjacent = filter newtest . filter (inBounds arr) $ adjDiagonal x
	in 
		if null adjacent 
			then let stack' = dropWhile (not . newtest) stack in --make sure there isn't an outdated element at the front
				if null stack'
					then boolArr
					else floodFillHelper test (head stack') arr  ( boolArr  V.// [(head stack',True)] ) (tail stack')
			else floodFillHelper test (head adjacent) arr ( boolArr  V.// [(head adjacent,True)] ) ( (tail adjacent) ++ stack)
--

adjCardinal :: Integral a => (a,a) -> [(a,a)]
adjCardinal (u,v) = [(u + 1,v),(u - 1,v),(u,v + 1),(u, v - 1)]
--

adjDiagonal :: Integral a => (a,a) -> [(a,a)]
adjDiagonal (u,v) = [(u + 1,v),(u - 1,v),(u,v + 1),(u, v - 1),(u + 1,v + 1),(u - 1,v - 1),(u - 1,v + 1),(u + 1, v - 1)]

--we are going to use this for oceans
test9 :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> V.Array (a,a) b
test9 = arrayMap (\x -> if x then 1 else 0) .  floodFill ( == 0) (0,0)

inBounds :: (V.Ix a,Ord a) =>  V.Array (a,a) b -> (a,a) -> Bool
inBounds arr (x,y)  = and [n_ <= x , x <= n' , m_ <= y , y <= m']
	where ((n_,m_),(n',m')) = V.bounds arr
--

findMax :: (V.Ix a,Enum a,Ord b) => V.Array (a,a) b -> ((a,a),b)
findMax arr = foldr (\n b -> foldr (\m c -> let d = ((n,m),(arr V.! (n,m))) in if (snd d) > (snd c) then d else c  ) b [m_..m'] ) ((n_,m_),(arr V.! (n_,m_))) [n_..n']
	where ((n_,m_),(n',m')) = V.bounds arr
--
--forest
test10 :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> V.Array (a,a) b
test10 = arrayMap (\x -> if ((x > 3) && (x < 13)) then 1 else 0)

--(r,g,b) --we assume each array has the same dimension
arrayToPPM :: (V.Ix b,Enum b,Num b,Show b,Show a,Num a,Ord a) => (V.Array (b,b) a,V.Array (b,b) a,V.Array (b,b) a) -> String
arrayToPPM (arr1,arr2,arr3) = "P3\n" ++ dim ++ "\n" ++ height ++ "\n" ++ arrList
	where
		((n_,m_),(n',m')) = V.bounds arr1
		dim = show (n' - n_ + 1) ++ " " ++ show (m' - m_ + 1)
		height = show $ maximum [height1,height2,height3]
		height1 = maximum . concat . arrayToList $ arr1
		height2 = maximum . concat . arrayToList $ arr2
		height3 = maximum . concat . arrayToList $ arr3
		arrList = unlines . map unwords . arrayToList . arrayMap show' . zip3Array $ (arr1,arr2,arr3)
		show' (x,y,z) = show x ++ " " ++ show y ++ " " ++ show z

zip3Array :: (V.Ix b,Enum b) => (V.Array (b,b) a,V.Array (b,b) a,V.Array (b,b) a) -> V.Array (b,b) (a,a,a)
zip3Array (arr1,arr2,arr3) = V.array ((n_,m_),(n',m')) [((x,y), ( arr1 V.! (x,y) , arr2 V.! (x,y) , arr3 V.! (x,y) )  ) | x <- [n_..n'] , y <- [m_..m'] ]
	where
		((n_,m_),(n',m')) = V.bounds arr1

greyToRGB :: V.Array (b,b) a -> (V.Array (b,b) a,V.Array (b,b) a,V.Array (b,b) a)
greyToRGB arr = (arr,arr,arr)

--we assume that alpha channels are between 0 and 1 -- 0 totally transparent and 1 totally opaque
--we are adding arr2 on top of arr1
addArrayAlpha ::  (V.Ix a,Enum a,Integral b) => Double ->  V.Array (a,a) b -> V.Array (a,a) b  -> V.Array (a,a) b
addArrayAlpha alpha2 arr2 arr1 = addArrays ( arrayMap (\x -> round $ (fromIntegral x) * alpha2) arr2 ) ( arrayMap (\x -> round $ (fromIntegral x) * (1 - alpha2)) arr1 )

ocean :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> ( V.Array (a,a) b , V.Array (a,a) b , V.Array (a,a) b)
ocean arr = ( constArray bnds 0 , constArray bnds 0 , arrayMap (\x -> if x then 15 else 0) . floodFill ( <= 1) (0,0) $ arr  )
	where
		bnds = V.bounds arr
--
forest :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> ( V.Array (a,a) b , V.Array (a,a) b , V.Array (a,a) b)
forest arr = ( constArray bnds 0 ,  arrayMap (\x -> if ((x > 3) && (x < 12)) then 15 else 0)  arr , constArray bnds 0  )
	where
		bnds = V.bounds arr
--
sand :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> ( V.Array (a,a) b , V.Array (a,a) b , V.Array (a,a) b)
sand arr = ( arrayMap (\x -> if ((x > 1) && (x <= 3)) then 15 else 0)  arr ,  arrayMap (\x -> if ((x > 1) && (x <= 3)) then 15 else 0)  arr ,arrayMap (\x -> if ((x > 1) && (x <= 3)) then 10 else 0)  arr  )
	where
		bnds = V.bounds arr

boolArrayNot :: (V.Ix a,Enum a) => V.Array (a,a) Bool -> V.Array (a,a) Bool
boolArrayNot = arrayMap not

boolArrayUnion :: (V.Ix a,Enum a) => V.Array (a,a) Bool -> V.Array (a,a) Bool -> V.Array (a,a) Bool
boolArrayUnion arr1 arr2 = V.array ((n_,m_),(n',m'))[((x,y), ( arr1 V.! (x,y) ) || ( arr2 V.! (x,y) ) ) | x <- [n_..n'] , y <- [m_..m'] ]
	where
		((n_,m_),(n',m')) = V.bounds arr1

boolArrayIntersection :: (V.Ix a,Enum a) => V.Array (a,a) Bool -> V.Array (a,a) Bool -> V.Array (a,a) Bool
boolArrayIntersection arr1 arr2 = V.array ((n_,m_),(n',m'))[((x,y), ( arr1 V.! (x,y) ) && ( arr2 V.! (x,y) ) ) | x <- [n_..n'] , y <- [m_..m'] ]
	where
		((n_,m_),(n',m')) = V.bounds arr1

lakes :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> ( V.Array (a,a) b , V.Array (a,a) b , V.Array (a,a) b)
lakes arr = (constArray bnds 0 , arrayMap (\x -> if x then 10 else 0) lks   , arrayMap (\x -> if x then 15 else 0) lks )
	where 
		bnds = V.bounds arr
		ocn = floodFill ( <= 1) (0,0) $ arr
		water = arrayMap (\x -> x <= 1) arr
		lks = boolArrayIntersection water (boolArrayNot ocn)

heightToColour :: (V.Ix a,Enum a,Integral a,Integral b) => V.Array (a,a) b -> ( V.Array (a,a) b , V.Array (a,a) b , V.Array (a,a) b)
heightToColour arr = dot3 (addArrayAlpha 0.5) ( lakes arr) . dot3 (addArrayAlpha 0.5) ( sand arr) . dot3 (addArrayAlpha 0.5 ) ( ocean arr) . dot3 (addArrayAlpha 0.5) (forest arr) $ (greyToRGB arr)


\end{code}

