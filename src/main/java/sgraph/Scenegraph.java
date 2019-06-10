package sgraph;


import org.joml.Matrix4f;
import org.joml.Vector3f;
import org.joml.Vector4f;
import util.*;
import raytracer.*;

import java.awt.*;
import java.io.IOException;
import java.util.*;
import java.util.List;

import static java.lang.Float.NaN;

/**
 * A specific implementation of this scene graph. This implementation is still
 * independent of the rendering technology (i.e. OpenGL)
 *
 * @author Amit Shesh
 */
public class Scenegraph<VertexType extends IVertexData> implements IScenegraph<VertexType> {
  /**
   * The root of the scene graph tree
   */
  protected INode root;

  /**
   * A map to store the (name,mesh) pairs. A map is chosen for efficient search
   */
  protected Map<String, util.PolygonMesh<VertexType>> meshes;

  /**
   * A map to store the (name,node) pairs. A map is chosen for efficient search
   */
  protected Map<String, INode> nodes;

  protected Map<String, String> textures;

  /**
   * The associated renderer for this scene graph. This must be set before
   * attempting to render the scene graph
   */
  protected IScenegraphRenderer renderer;

  //textureimage for raytracing
  private Map<String, TextureImage> textureImages;


  public Scenegraph() {
    Matrix4f mv;
    root = null;
    meshes = new HashMap<String, util.PolygonMesh<VertexType>>();
    nodes = new HashMap<String, INode>();
    textures = new HashMap<String, String>();
    textureImages = new HashMap<String, TextureImage>();
    try {
      textureImages.put("white", new TextureImage("textures/white.png", "png", "white"));
    }catch (IOException e){
      e.printStackTrace();
    }
  }

  public void dispose() {
    renderer.dispose();
  }

  /**
   * Sets the renderer, and then adds all the meshes to the renderer. This
   * function must be called when the scene graph is complete, otherwise not all
   * of its meshes will be known to the renderer
   *
   * @param renderer The {@link IScenegraphRenderer} object that will act as its
   *                 renderer
   */
  @Override
  public void setRenderer(IScenegraphRenderer renderer) throws Exception {
    this.renderer = renderer;

    //now add all the meshes
    for (String meshName : meshes.keySet()) {
      this.renderer.addMesh(meshName, meshes.get(meshName));
    }

    //pass all the texture objects
    for (Map.Entry<String, String> entry : textures.entrySet()) {
      renderer.addTexture(entry.getKey(), entry.getValue());
    }

  }


  /**
   * Set the root of the scenegraph, and then pass a reference to this scene
   * graph object to all its node. This will enable any node to call functions
   * of its associated scene graph
   */

  @Override
  public void makeScenegraph(INode root) {
    this.root = root;
    this.root.setScenegraph(this);

  }

  /**
   * Draw this scene graph. It delegates this operation to the renderer
   */
  @Override
  public void draw(Stack<Matrix4f> modelView) {
    if ((root != null) && (renderer != null)) {
      List<Light> listOfLights = root.getLightsInView(modelView);
      renderer.initLightsInShader(listOfLights);
      renderer.draw(root, modelView);
    }
  }


  @Override
  public void addPolygonMesh(String name, util.PolygonMesh<VertexType> mesh) {
    meshes.put(name, mesh);
  }


  @Override
  public void animate(float time) {

  }

  @Override
  public void addNode(String name, INode node) {
    nodes.put(name, node);
  }


  @Override
  public INode getRoot() {
    return root;
  }

  @Override
  public Map<String, PolygonMesh<VertexType>> getPolygonMeshes() {
    Map<String, util.PolygonMesh<VertexType>> meshes = new HashMap<String, PolygonMesh<VertexType>>(this.meshes);
    return meshes;
  }

  @Override
  public Map<String, INode> getNodes() {
    Map<String, INode> nodes = new TreeMap<String, INode>();
    nodes.putAll(this.nodes);
    return nodes;
  }

  @Override
  public void addTexture(String name, String path) {
    textures.put(name, path);
  }

  /**
   * Compute the ray starting from the eye and passing through this pixel in view coordinates.
   * Pass the ray to a function called “raycast” that returns information about ray casting.
   * Write the color to the appropriate place in the array.
   * @param w
   * @param h
   * @param modelView
   * @return
   */
  @Override
  public Color[][] raytrace(int w, int h, Stack<Matrix4f> modelView) {
    Color[][] colors = new Color[w][h];
    Vector4f eye = new Vector4f(0, 0, 0, 1);
    float fieldOfView = 120.0f;

    for (int i = 0; i < w; i++) {
      for (int j = 0; j < h; j++) {
        // pixel in view coordinates
        float deltaX = i - w / 2 - eye.x;
        float deltaY = j - h / 2 - eye.y;
        float deltaZ = (float) (-0.5f * h / Math.tan((float) Math.toRadians(fieldOfView) / 2));

        Ray ray = new Ray(eye, new Vector4f(deltaX, deltaY, deltaZ, 0));
        Vector4f bgColor = raycast(ray, modelView,0,true);
        colors[i][j] = new Color(bgColor.x, bgColor.y, bgColor.z);
      }
    }

    return colors;
  }

  /**
   * 1. Pass on the ray and the modelview stack to the root of the scenegraph.
   * 2. If the ray did hit anything, call a function shade and pass to it all the relevant information from the hit record that the root returns, that will be useful in calculating the color of this pixel.
   * 3. If the ray did not hit anything, return the background color.
   * @param ray incoming ray
   * @param modelView view transform
   * @return color
   */
  @Override
  public Vector4f raycast(Ray ray, Stack<Matrix4f> modelView,int count,boolean isRayInAir) {
    Vector4f bgColor = new Vector4f(0, 0, 0, 1);
    if(count<5){
      HitRecord hitRecord = getRoot().calcIntersect(ray, modelView);
      if (hitRecord.isHit()) {
        List<Light> lights = this.getRoot().getLightsInView(modelView);
        Vector4f color = this.shade(hitRecord, lights,ray,count,modelView,isRayInAir);
        bgColor.set(color.x, color.y, color.z, 0);
      }
    }
    return bgColor;
  }

  /**
   *
   *
   * @param ray incoming ray
   * @param modelView view transform
   * @param count ray bounce count
   * @param reflection reflection index
   * @param rcolor reflection color
   * @return
   */
  private Vector4f rayCastForReflection(Ray ray, Stack<Matrix4f> modelView,int count,float reflection,Vector4f rcolor){

      rcolor = raycast(ray,modelView,count,false);
      rcolor = rcolor.mul(reflection);
      return rcolor;
  }

  /**
   *
   * @param ray incoming ray
   * @param modelView view transform
   * @param count bounce count
   * @param refraction transparency index index
   * @param tcolor  transparency color
   * @param isRayInAir ray in air or not
   * @return color
   */
  private Vector4f rayCastForRefraction(Ray ray, Stack<Matrix4f> modelView,int count,float refraction,Vector4f tcolor,boolean isRayInAir){

    tcolor = raycast(ray,modelView,count,isRayInAir);
    tcolor = tcolor.mul(refraction);
    return tcolor;
  }


  /**
   *  shade pixel color
   * @param hitRecord
   * @param lights all lights in scene
   * @return color
   */

//  vec3 lightVec,viewVec,reflectVec;
//  vec3 normalView;
//  vec3 ambient,diffuse,specular;
//  float nDotL,rDotV;
//
//
//  fColor = vec4(0,0,0,1);
//
//    for (int i=0;i<numLights;i++)
//  {
//    if (light[i].position.w!=0)
//      lightVec = normalize(light[i].position.xyz - fPosition.xyz);
//    else
//      lightVec = normalize(-light[i].position.xyz);
//
//    vec3 spotdirection = normalize(light[i].spotdirection.xyz);
//    if (dot(-lightVec,spotdirection)<light[i].cosSpotCutoff)
//      continue;
//
//    vec3 tNormal = fNormal;
//    normalView = normalize(tNormal.xyz);
//
//
//    nDotL = dot(normalView,lightVec);
//
//    viewVec = -fPosition.xyz;
//    viewVec = normalize(viewVec);
//
//    reflectVec = reflect(-lightVec,normalView);
//    reflectVec = normalize(reflectVec);
//
//    rDotV = max(dot(reflectVec,viewVec),0.0);
//
//    ambient = material.ambient * light[i].ambient;
//    diffuse = material.diffuse * light[i].diffuse * max(nDotL,0);
//    if (nDotL>0)
//      specular = material.specular * light[i].specular * pow(rDotV,material.shininess);
//    else
//      specular = vec3(0,0,0);
//    fColor = fColor + vec4(ambient+diffuse+specular,1.0);
//  }
//  fColor = fColor * texture(image,fTexCoord.st);
  public Vector4f shade(HitRecord hitRecord, List<Light> lights, Ray currentRay, int bounceCount, Stack<Matrix4f> mv,boolean isRayInAir) {
    Vector4f fcolor = new Vector4f(0, 0, 0, 0);
    Vector4f intersectPosition = new Vector4f(hitRecord.getIntersectionP_view());
    Vector3f normalView = new Vector3f(getXYZ(hitRecord.getNormal()));
    normalView = normalView.normalize();
    Material material = hitRecord.getMaterial();
    float a = material.getAbsorption();
    float r = material.getReflection();
    float t = material.getTransparency();
    //for each light

    for (int i = 0; i < lights.size(); i++) {
      //get lightvec
      Light light = lights.get(i);
      Vector3f lightVec;
      if (light.getPosition().w!=0) {
        lightVec = new Vector3f(getXYZ(light.getPosition())).sub(getXYZ(intersectPosition));
      }
      else {
        lightVec = new Vector3f(getXYZ(light.getPosition())).negate();
      }
      lightVec = lightVec.normalize();
      Vector3f lightXYZ = new Vector3f(lightVec);

      //System.out.println(normalView);
      Vector3f lightVecNegate = new Vector3f();
      lightVec.negate(lightVecNegate);
      Vector3f spotDirection = new Vector3f(getXYZ(light.getSpotDirection())).normalize();
      float cosAngle = (float) Math.cos(Math.toRadians(light.getSpotCutoff()));
     // System.out.println(spotDirection.dot(lightVecNegate));
      if (spotDirection.dot(lightVecNegate)<cosAngle)
        continue;

      float nDotL = normalView.dot(lightVec);

      /**
       * Shadow
       */
      boolean isShadow = false;
      Vector4f shadowRayDir = new Vector4f(new Vector3f(lightXYZ), 0);
      Vector4f fudge = new Vector4f(shadowRayDir).mul(0.01f);
      Vector4f startPoint = new Vector4f(intersectPosition);
      startPoint.add(fudge);

      Ray ray = new Ray(startPoint, shadowRayDir);
      HitRecord _hitRecord = getRoot().calcIntersect(ray, mv);
      if (_hitRecord.isHit()) {
        if(light.getPosition().w == 1){
          float disToLight = startPoint.distance(light.getPosition());
          float disToIntersection = startPoint.distance(_hitRecord.getIntersectionP_view());
          if(disToIntersection<=disToLight)
            isShadow = true;
        }else{
          isShadow = true;
        }

      }
//      if (light.getPosition().w != 0) {
//        if (_hitRecord.isHit()) {
//          isShadow = true;
//        }
//      } else {
//        float dist = light.getPosition().distance(startPoint);
//        if (_hitRecord.isHit() && dist - _hitRecord.getT() > 0.001f) {
//          isShadow = true;
//        }
//      }

  //    System.out.println("Shadow: " + isShadow);


      Vector3f viewVec = new Vector3f(getXYZ(intersectPosition)).negate();
      viewVec = viewVec.normalize();

      Vector3f reflectVec = new Vector3f(lightVec.negate().reflect(normalView));
      reflectVec = reflectVec.normalize();

      double rDotV = Math.max(reflectVec.dot(viewVec),0.0);
      //ambient

      Vector3f materialAmbient = getXYZ(material.getAmbient());
      Vector3f lightAmbient = new Vector3f(light.getAmbient());
      Vector3f ambient = new Vector3f(materialAmbient);
      if(!isShadow && isRayInAir)
        ambient.mul(new Vector3f(lightAmbient));
      //diffuse
      Vector3f materialDiffuse = getXYZ(material.getDiffuse());
      Vector3f lightDiffuse = new Vector3f(light.getDiffuse());
      Vector3f diffuse = new Vector3f(materialDiffuse).mul(new Vector3f(lightDiffuse))
                      .mul(Math.max(0, nDotL));


      //specular
      Vector3f specular;
      if (nDotL > 0) {
        Vector3f materialSpecular = getXYZ(material.getSpecular());
        Vector3f lightSpecular = new Vector3f(light.getSpecular());
        specular = new Vector3f(materialSpecular).mul(new Vector3f(lightSpecular));
        specular = specular.mul((float) Math.pow(rDotV, material.getShininess()));
      }
      else {
        specular = new Vector3f(0,0,0);
      }
      if(isShadow || !isRayInAir){
        diffuse = new Vector3f(0,0,0);
        specular = new Vector3f(0,0,0);
      }
      //add ambient, diffuse, specular
      Vector3f addup = new Vector3f().add(ambient).add(diffuse).add(specular);
      fcolor = fcolor.add(new Vector4f(addup,1));

    }

    //add texture
    applyTexture(fcolor,hitRecord);

    //reflection//////////////////////////////////////////////////////////
    Vector4f rcolor = new Vector4f(0, 0, 0, 0);
    if(r!=0 && isRayInAir){
      Vector3f incomingRayDir = new Vector3f(currentRay.getDirection().x,currentRay.getDirection().y,currentRay.getDirection().z);
      Vector3f reflectRayDir = new Vector3f();
      incomingRayDir.reflect(normalView,reflectRayDir);
      Vector4f reflectStartPoint = new Vector4f(intersectPosition);
      reflectStartPoint.add(new Vector4f(reflectRayDir.x,reflectRayDir.y,reflectRayDir.z,0).mul(0.01f));
      Ray reflectRay = new Ray(reflectStartPoint,new Vector4f(reflectRayDir.x,reflectRayDir.y,reflectRayDir.z,0));
      rcolor = rayCastForReflection(reflectRay,mv,bounceCount+1,r,rcolor);

    }
    //reflection index multiplied in raycast
    fcolor.mul(a).add(rcolor);


    //refraction/////////////////////////////////////////////////////////////
    Vector4f tcolor = new Vector4f(0,0,0,1);
    if(t!=0){
      Vector3f refractionRayDir = new Vector3f();
      Vector3f incomingRayDir = getXYZ(currentRay.getDirection());
      Vector3f intersectionNormal = new Vector3f(normalView);
      Vector4f intersectionPoint = new Vector4f(intersectPosition);
      float cosThetaI = -1f*intersectionNormal.dot(incomingRayDir);
      float sinThetaISqr = 1f - cosThetaI*cosThetaI;
      float materialIndex = 1.05f;
      float ratio = 0;
      if (isRayInAir) {
        // ray hit in air
        ratio = 1f / materialIndex;
      } else {
        // ray hit in object
        ratio = materialIndex / 1f;
        intersectionNormal.negate();
      }
      float sinThetaRSqr = ratio * ratio * sinThetaISqr;
      float cosThetaRSqr = 1-sinThetaRSqr;

      // ray inside object can go out
      if (cosThetaRSqr >= 0f) {
        float cosThetaR = (float) Math.sqrt(cosThetaRSqr);
        refractionRayDir = (new Vector3f(intersectionNormal).mul(cosThetaI).add(incomingRayDir)).mul(ratio)
                .sub(new Vector3f(intersectionNormal).mul(cosThetaR));

        Vector4f refractionDir4f = new Vector4f(refractionRayDir, 0).normalize();
        Ray refractionRay = new Ray(intersectionPoint.add(new Vector4f(refractionDir4f).mul(0.01f)),
                refractionDir4f);
        tcolor = rayCastForRefraction(refractionRay, mv, bounceCount+1,t,tcolor,!isRayInAir);
      }
      //refraction index multiplied in raycast
      fcolor.add(tcolor);
    }
    normalizeColor(fcolor);
    return fcolor;
  }

  /**
   * get x,y,z of a vector4f
   * @param v vector
   * @return vector3f that has same x,y,z
   */
  private Vector3f getXYZ(Vector4f v){
    return new Vector3f(v.x,v.y,v.z);
  }

  /**
   *  apply texture color to fcolor
   * @param fcolor color without texture
   * @param record hitrecord
   */
  public void applyTexture(Vector4f fcolor, HitRecord record){
    Vector4f imageColor;
    if(record.getTextureName().equals("")){
      imageColor = textureImages.get("white").getColor(record.getTextureX(), record.getTextureY());
    }else{
      imageColor = textureImages.get(record.getTextureName()).getColor(record.getTextureX(), record.getTextureY());
    }

    fcolor = fcolor.mul(imageColor);
  }

  /**
   *
   * @param name name of image
   * @param path path to image
   *             set image map
   */
  @Override
  public void setTextureImages(String name, String path){
    String imageFormat = path.substring(path.indexOf('.') + 1);
    try{
      TextureImage tximage = new TextureImage(path,imageFormat,name);
      this.textureImages.put(name,tximage);
    }catch (IOException e){
      e.printStackTrace();
    }

  }

  /**
   *
   * @param color
   * normalize rgb value of color
   */
  private void normalizeColor(Vector4f color) {
    if (color.x > 1) {
      color.x = 1;
    }
    if (color.x < 0) {
      color.x = 0;
    }
    if (color.y > 1) {
      color.y = 1;
    }
    if (color.y < 0) {
      color.y = 0;
    }
    if (color.z > 1) {
      color.z = 1;
    }
    if (color.z < 0) {
      color.z = 0;
    }
  }
}
