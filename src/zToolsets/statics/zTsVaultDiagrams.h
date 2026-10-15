#pragma once
#include <zspace/zToolsets/statics/zTsVault.h>
namespace zSpace {
struct zVaultDiagramData {
    zObjectMesh form, force;
    zIntArray formToForce, forceToForm;
    zVaultDiagramSettings settings;
    zVaultDiagramReport report;
};
zVaultDiagramData buildVaultDiagrams(zObjectMesh& result, const zIntArray& supports,
    const zIntArray& edges, const zDoubleArray& q, const zVaultDiagramSettings& settings);
}
