/* variant.c -- routines to convert other variant formats to a generic
 *              variant structure */

#include "common.h"
#include "variant.h"

struct allele  *alleleClip(struct allele *allele, int sx, int ex)
/* clip allele to be inside region defined by sx..ex.  Returns 
 * pointer to new allele which should be freed by alleleFree, or variantFree
 */
{
struct variant *oldVariant = allele->variant;
int start = oldVariant->chromStart;
int end = oldVariant->chromEnd;
int oldVariantWidth = end - start;
int delFront = 0;
int delRear = 0;

if (start < sx)
    {
    if (oldVariantWidth != allele->length)	 /* FIXME */
	errAbort("cannot clip alleles that are a different length than variant region");
    delFront = sx - start;
    start = sx;
    }

if (end > ex)
    {
    if (oldVariantWidth != allele->length)	 /* FIXME */
	errAbort("cannot clip alleles that are a different length than variant region");
    delRear = end - ex;
    end = ex;
    }

struct variant *newVariant;
AllocVar(newVariant);
newVariant->chrom = cloneString(oldVariant->chrom);
newVariant->chromStart = start;
newVariant->chromEnd = end;
newVariant->numAlleles = 1;
struct allele *newAllele;
AllocVar(newAllele);
newVariant->alleles = newAllele;
newAllele->variant = newVariant;
newAllele->length = allele->length - delRear - delFront;
assert(newAllele->length > 0);
newAllele->sequence = cloneString(&allele->sequence[delFront]);
newAllele->sequence[newAllele->length] = 0;   // cut off delRear part

return newAllele;
}

static char *makeDashes(int count)
{
char *ret = needMem(count + 1);
char *ptr = ret;

while(count--)
    *ptr++ = '-';

return ret;
}

struct variant *variantFromPgSnp(struct pgSnp *pgSnp)
/* convert pgSnp record to variant record */
{
struct variant *variant;

// this is probably the wrong way to do this.  Alleles in
// variant should be their size in query bases
int alleleLength = pgSnp->chromEnd - pgSnp->chromStart;

// We have a new variant!
AllocVar(variant);
variant->chrom = cloneString(pgSnp->chrom);
variant->chromStart = pgSnp->chromStart;
variant->chromEnd = pgSnp->chromEnd;
variant->numAlleles = pgSnp->alleleCount;

// get the alleles.
char *nextAlleleString = pgSnp->name;
int alleleNumber = 0;
for( ; alleleNumber < pgSnp->alleleCount; alleleNumber++)
    {
    if (nextAlleleString == NULL)
	errAbort("number of alleles in pgSnp doesn't match number in name");
    
    char *thisAlleleString = nextAlleleString;

    // advance pointer to next variant string
    // probably there's some kent routine to do this behind the curtain
    nextAlleleString = strchr(thisAlleleString, '/');
    if (nextAlleleString)	 // null out '/' and move to next char
	{
	*nextAlleleString = 0;
	nextAlleleString++;
	}

    // this check probably not right, could be different per allele
    int alleleStringLength = strlen(thisAlleleString);
    if (alleleStringLength != alleleLength)
	{
	// check for special case of single '-'
	if (sameString("-", thisAlleleString))
	    {
	    thisAlleleString = makeDashes(alleleLength);
	    alleleStringLength = alleleLength;
	    }
	else
	    errAbort("length of allele number %d is %d, should be %d", 
		alleleNumber, alleleStringLength, alleleLength);
	}

    // we have a new allele!
    struct allele *allele;
    AllocVar(allele);
    slAddHead(&variant->alleles, allele);
    allele->variant = variant;
    allele->length = alleleStringLength; 
    toLowerN(thisAlleleString, alleleStringLength);
    allele->sequence = cloneString(thisAlleleString);
    }

slReverse(&variant->alleles);

return variant;
}

