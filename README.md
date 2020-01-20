# MorphTrans: word-word translation in morphologically complex languages

MorphTrans is a project dedicated to investigating the machine translation subtask word-to-word translation.
If all languages were as uninteresting morphologically as English and Mandarin, there would be no reason to spend time on this subtask.
But for the thousands of languages in the world with complex morphology, including many low-resource languages, word-to-word translation has much of the complexity of sentence-to-sentence translation.

In this project, we focus on the subtask within this subtask of translation between morphologically complex languages that are closely related, specifically on the roughly 15 languages within the South Semitic group spoken in Ethiopia and Eritrea, usually referred to as "Ethiopian Semitic" or "Ethio-Semitic" (somewhat confusingly because they are spoken in two countries; in what follows, we'll refer to them as "ES languages"). Languages we are working with include Amharic (Am, አማርኛ), Tigrinya (Tn, ትግርኛ), and Kistane (Ks, ክስታንኘ, Kistaninya, Soddo).

## Data files
Translation data files are in subdirectories in Data/, one for each
language pair.
Current folders include:

* AmKs
    * Amharic and Kistane
    * Kistane roots and their Amharic translations were extracted from the Kistane-Amharic-English dictionary compiled by the Kistane People's Development Association and digitized by Ermiyas Weldeyohanes.
    * Word forms were generated by [HornMorpho](https://github.com/hltdi/HornMorpho).

Translation files end with the extension .tr.

### Translation file format
Each line in a translation data (.tr) file has the following form
```
source_word;target_words;;source_root_asp;target_root_asp;gramfeats
```
where `target_words` is one or more words in the target language (separated by commas), `source_root_asp` and `target_root_asp` are combinations of roots and aspect/voice categories (see below for details for ES languages), and `gramfeats` is a comma-separated list of grammatical feature abbreviations.

### Romanization
ES languages are written using the syllabic Ge'ez script but are romanized for this project because morphological generalizations are more easily made at the level of the individual segment.

We use the ASCII romanization scheme that is also used in [HornMorpho](https://github.com/hltdi/HornMorpho).

Important: there are some two-character segments (especially those with 'W' as the second character), so ``list(word)`` won't work for extracting the segments in a word string. Use the function ``get_words()`` in ``Data/proc.py`` instead; it segments word strings correctly.

### Aspect/voice categories
In ES languages, each verb root, consisting of a sequence of consonants, can appear in up to ten different patterns, which can be viewed as combinations of various aspect and voice features.
Note that these are **formal** categories; for some roots the labels may not correspond closely to their function (see the next paragraph).
For example, from the root `<sbr>` we can generate the "passive-iterative" form (abbreviated "psit"), which seems to be common to all of the languages (Am *tesebab\_ere*, Tn *tesebabere*, Ks *tesbab\_ero*).

Between languages, rather than translating roots, it is more reasonable to translate combinations of roots and aspect/voice categories, in effect treating these combinations as the lexemes.
This is because the actual function of the categories varies across roots and between languages.
For example, the Am root `<drg>` 'do' fails to appear in the basic aspect/voice category; its basic form is instead the "transitive" form *ader_ege*. But this translates in Kistane to an unrelated root, `<qnh.B>`, in the basic category, *qin_aw*.

### Grammatical features

Each combination of root and aspect/voice category in the source language is combined with multiple sets of grammatical features to generate the final word forms in a translation file.
Some of the features apply to all words, while those related to subject and object person, number, and gender depend on the valency of the combination of root and aspect/voice category.

Since valency is not provided in the available dictionaries, we need to extract it automatically. For Amharic and Tigrinya, we have the useless lists of unique word types gathered by Biniam Gebremichael: https://www.cs.ru.nl/~biniam/geez/crawl.php. Using these words and the morphological analyzer of [HornMorpho](https://github.com/hltdi/HornMorpho), it is possible to classify many root/aspect-voice combinations by their valency. For Amharic we have the following categories:

* Subject-less (only 3 masculine singular subjects, no objects)
* Impersonal (only 3 masculine singular subjects, any objects)
* Intransitive (any subjects, no objects)
* Transitive (any subjects, any objects)

### Specific files
* 1.py
    * 23,908 translation pairs
    * 712 root/aspect-voice categories
    * Gemination is indicated for both Amharic and Kistane words (represented by '_' following the geminated consonant
    * Epenthetic (sixth-order) vowels are not indicated
    * Grammatical features
        * perfective vs. imperfective
        * main vs. relative/subordinate
        * affirmative vs. negative
        * subjects (depending on valency): 3 singular masculine, 3 singular feminine, 1 plural
        * objects (depending on valency): 1 singular, 3 singular masculine, 2 plural
    * Missing:
        * Jussive/imperative
        * Applicative objects (ፈረደለት fer_edel_et)
        * Ks "impersonal", usually corresponding to Am passive (ይብሉት yblut = Am ይባላል yb_alal)
    




