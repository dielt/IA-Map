\begin{code}


module Util.Base where


import Data.Numbers.Primes
import Data.Maybe
import Data.Monoid
import Data.Tree
import Data.List
import Control.Monad
import Control.Monad.Trans.State
import qualified Control.Category as C
import Control.Monad
import Control.Applicative

\end{code}

We can alway split this up into seperate files
e.g. IAStateT IACircuit or similar

Basic utility stuff,

\begin{code}
infixr 8 .:
(.:) :: (c->d) ->(a->b->c) -> (a->b->d)
(.:) = (.) . (.)

infixr 8 .::
(.::) :: (d->e) -> (a->b->c->d) -> (a->b->c->e)
(.::) = (.) . (.) . (.)

fst3 (x,_,_) = x
snd3 (_,x,_) = x
thd3 (_,_,x) = x

appFst f (x,y) = (f x ,y)
appSnd f (x,y) = (x, f y)

if' :: Bool -> a -> a -> a
if' True  x _ = x
if' False _ y = y

--sort of like guard or when but with values
ifM :: MonadPlus m => Bool -> m a -> m a
ifM True  x = x
ifM False _ = mzero 

notNothing :: Maybe a -> Bool
notNothing = not . isNothing

uncurry3 :: (a -> b -> c -> d) -> (a,b,c) -> d
uncurry3 f (a,b,c) = f a b c

uncurry4 :: (a -> b -> c -> d -> e) -> (a,b,c,d) -> e
uncurry4 f (a,b,c,d) = f a b c d


--I don't think the type system will allow for a generalized version
eat1Arg f = \a -> f
eat2Arg f = \a b -> f
eat3Arg f = \a b c -> f
eat4Arg f = \a b c d -> f
eat5Arg f = \a b c d e -> f
eat6Arg f = \a b c d e g -> f


head' = listToMaybe
tail' xs = if null xs then Nothing else Just $ tail xs

emptyTail [] = []
emtpyTail xs = tail xs

appHead _ [] = []
appHead f (x:xs) = (f x) : xs

appTail _ [] = []
appTail f (x:xs) = x : (map f xs)

--again I assume that there is a prelude function for this
deleteAll :: Eq a => a -> [a] -> [a]
deleteAll = filter . (/=)



\end{code}


Some more specialized functions for dealing with lists

\begin{code}

listFstFilter :: [(Bool,a)] -> [a]
listFstFilter xs = foldr (\x list -> if fst x then snd x : list else list ) [] xs

--based on the last of the tuple
findLowSnd2 :: Ord b => [(a,b)] -> Maybe (a,b)
findLowSnd2 [] = Nothing
findLowSnd2 (x:xs) = Just $ foldr (\a b -> if (snd a) <= (snd b) then a else b ) x xs
 
findHighSnd2 :: Ord b => [(a,b)] -> Maybe (a,b)
findHighSnd2 [] = Nothing
findHighSnd2 (x:xs) = Just $ foldr (\a b-> if snd a > snd b then a else b) x xs

findLowThd3 :: Ord c => [(a,b,c)] -> Maybe (a,b,c)
findLowThd3 [] = Nothing
findLowThd3 (x:xs) = Just $ foldr (\a b -> if (thd3 a) <= (thd3 b) then a else b ) x xs


checkFst3 :: Eq a => a -> [(a,b,c)] -> Bool
checkFst3 a list = or $ map ((a ==) . fst3) list
checkSnd3 :: Eq b => b -> [(a,b,c)] -> Bool
checkSnd3 a list = or $ map ((a ==) . snd3) list
checkThd3 :: Eq c => c -> [(a,b,c)] -> Bool
checkThd3 a list = or $ map ((a ==) . thd3) list

\end{code}


\begin{code}

listOr a b = if null a then b else a

emptyList :: a -> [a]
emptyList a = tail [a]

\end{code}



\begin{code}


joinTuple :: Monad m => m (m a , m b) -> ( m a , m b ) 
joinTuple x = let
	a' = x >>= (\(a,b) -> a )
	b' = x >>= (\(a,b) -> b )
	in (a',b')

--

appTuple1 :: (a->b) -> (a,a) -> (b,b) 
appTuple1 f (x,y) = (f x , f y)

appTuple2 :: (a->b->c) -> (a,a) -> (b,b) -> (c,c)
appTuple2 f (x,y) (z,w) = (f x z , f y w)

addTuple = appTuple2 (+)

foldTuple :: (a -> b -> b) -> (b,b) -> [(a,a)] -> (b,b)
foldTuple f = foldr ( appTuple2 f )

sumTuple = foldTuple (+) (0,0)

\end{code}


More Mathematical style sums and products, sequence style
\begin{code}

sumSeq :: (Enum a,Num b) => (a -> b) -> a -> a -> b
sumSeq f i n = sum [(f x) | x <- [i..n] ] 

prodSeq :: (Enum a,Num b) => (a -> b) -> a -> a -> b
prodSeq f i n = product [(f x) | x <- [i..n] ] 

foldSeq :: (Enum a) => (b -> c -> c) -> c -> (a -> b) -> a -> a -> c
foldSeq g z f i n = foldr g z [f x | x <- [i..n]]


\end{code}






\begin{code}

--forces a function to act per line on streaming textual input
eachLine :: (String -> String) -> (String -> String)
eachLine f = unlines . map f . lines

\end{code}


This should go into some sort of data.base thing
\begin{code}

newtype MTrue = MTrue { getTrue :: Bool }

instance Monoid MTrue where
	mempty = MTrue $ True
	mappend x y = MTrue $ (getTrue x) && (getTrue y) 

newtype MFalse = MFalse { getFalse :: Bool }

instance Monoid MFalse where
	mempty = MFalse $ False
	mappend x y = MFalse $ (getFalse x) || (getFalse y)



\end{code}



Note these hash functions are not at all garrenteed to have any nice properties, in regards to the distribution of numbers
\begin{code}

simpleHash1 x = sum ( map (x^) (take ((x `mod` 13) + 1) primes)) + foldl' (\y z -> y + foldl' (\a b -> a + b^z ) 0 l ) 0 (take ((x `mod` 7) + 1) primes)
	where l = (map (\a -> if a < 0 then a^2 else a) [(x - 10)..(x)])

simpleHash :: Int -> Int -> Int
simpleHash x r = ( (randomVals !! (x^5 `mod` 200) )+( (randomVals !! (x^2 `mod` 200) ) + simpleHash1 ( (randomVals !! (x^3 `mod` 200) ) + x * ( (simpleHash1 x) `mod` (r + 7) )))) `mod` r

		
		
		
		
--from random.org
randomVals = [24004,-53750,91624,-75456,72735,55191,97125,54032,-50418,96909,75604,34169,-98369,269,26464,-11127,-43890,-87904,10147,84452,58307,70985,-58964,40626,-47841,-60653,5949,95884,-89205,19815,-41635,-73069,-35471,-36051,52406,30699,20418,-32563,-83836,21416,57529,40825,-50534,3053,-71408,-55690,-16883,-596,82553,-11125,59738,-41558,11359,20723,24285,-86367,-24626,-45773,-60290,-12904,-97202,-52044,-9235,-65867,83827,30741,93542,-56142,-37966,58907,-82358,-47826,-48406,-31531,-93437,-82564,-46223,-10604,-90944,-47064,54811,68207,11020,-93973,-93557,52385,42780,-86033,25373,-19220,-59791,-91075,-76813,44626,-22723,-31659,-69648,86060,44707,53976,-9637,27690,36483,-25734,3099,95872,-46543,-43458,68681,879,-11768,12147,58683,-38319,27814,-91117,98377,-57268,-41117,42845,29051,31608,-56337,79017,-31081,72587,-56112,-52905,42384,94039,63416,20948,21711,-82247,61641,-12686,-68746,36241,93852,19659,13935,-96584,16463,15717,-93809,-27916,58326,-60303,-63261,-62583,40729,-18692,53426,71373,-82576,37777,-88555,43836,-24176,-45470,36527,-36900,-66272,-25470,-56319,40135,76581,79033,-27724,72328,68589,90846,22765,-72927,19524,-95788,31714,8371,-99996,-84982,46310,-21959,8066,-71515,17544,41609,-81798,1725,-43131,43121,-36255,-24926,-33633,-42816,-81578,43575,57262,4526,-98573,66039]

simpleHashList :: [Int] -> Int -> Int
simpleHashList [] r = 777777777 `mod` r 
simpleHashList (x:xs) r = ( (randomVals !! ((x + 17) `mod` 200) ) + simpleHash x ( r + 7 ) + simpleHash (randomVals !! (x `mod` 200) ) r^3 + (simpleHashList xs (r + 3)^3)  ) `mod` r


\end{code}




















Merged stuff

\begin{code}

--stuff below should probably go into some sort of parser yModule

formatNames :: [String] -> String
formatNames [] = ""
formatNames (x:xs) = let
		counter stuff def = foldr (\a b -> if a == def then b+1 else b ) 1 stuff
		count =  counter xs x
		remaining = deleteAll x xs
	in
		(if count > 1
			then (show count) ++ " " ++ x ++ "s"
			else if (head x) `elem` upperCaseLetters
				then x
				else if (head x) `elem` vowels
					then "an " ++ x
					else "a " ++ x 
	) ++ (
		if null remaining
			then "" 
			else if null . deleteAll (head remaining) . tail $ remaining
				then ", and " 
				else ", "
		)  ++ (formatNames remaining)


upperCaseLetters :: [Char]
upperCaseLetters =
	"ABCDEFGHIJKLMNOPQRSTUVWXYZ"

vowels :: [Char]
vowels =
	"aeoiuAEOIU"
	
letters :: [Char] -- in alphebetic order
letters =
	"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

numbers :: [Char]
numbers = 
	"1234567890"

punctuation :: [Char]
punctuation = 
	"-_,.;:?"

characters = punctuation ++ numbers ++ letters


\end{code}




