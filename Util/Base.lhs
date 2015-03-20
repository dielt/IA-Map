\begin{code}


module Util.Base where


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




